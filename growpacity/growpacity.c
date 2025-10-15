#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

char *Array1D(int, int);
char ***Array3D(int, int, int, int);
void ReadOpacityData();
void ReadArray(char *, double **, int *, int);
double EvaluateRosselandOpacityArray(double q, double amax, double T);
double EvaluatePlanckOpacityArray(double q, double amax, double T);

static double *q_arr, *log_amax_arr, *log_T_arr;
static double ***log_kR_arr, ***log_kP_arr;
static int Nq, Namax, NT;
static double d_q_arr, dlog_amax_arr, dlog_T_arr;

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

char ***Array3D(int N1, int N2, int N3, int size) {
	// allocates a 3D array of N1 x N2 x N3 elements of given size (in bytes),
	// returns pointer to the array
	int i, j;
	char ***arr;
	arr = (char ***) malloc(N1 * sizeof(char **));
	for (i=0; i<N1; i++) {
		arr[i] = (char **) malloc(N2 * sizeof(char *));
		for (j=0; j<N2; j++) {
			arr[i][j] = (char *) malloc(N3 * size);
		}
	}
	return arr;
}

void ReadArray(char *flnm, double **arr, int *Npts, int use_log) {
	/*
	Reads a 1D array from file flnm.
	The first line of file should contain number of points.
	If use_log is set, takes log of values read.

	Allocates memory for array, sets Npts,
	returns pointer to array and number of points via arguments.
	*/

	FILE *f;

	f = fopen(flnm, "r");
	fscanf(f, "%d\n", Npts);

	*arr = (double *) malloc(*Npts * sizeof(double));

	int i;
	double tmp;
	for (i=0; i<*Npts; i++) {
		fscanf(f, "%le\n", &tmp);
		if (use_log) tmp = log(tmp);
		(*arr)[i] = tmp;
	}

	fclose(f);
}

void ReadOpacityData() {
	/*
	Reads opacity data from files. Also allocates memory for arrays.
	Files:
		q.dat       - first line number of points, then q values
		amax_um.dat - first line number of points, then amax values in microns
		T_K.dat     - first line number of points, then T values in K
		kR_cm2g.dbl - binary file with Rosseland mean opacities
		kP_cm2g.dbl - binary file with Planck mean opacities
	All arrays are in log space except q, indexing is (q, amax, T)
	*/

	// read grid arrays
	ReadArray("q.dat", &q_arr, &Nq, 0);
	ReadArray("amax_um.dat", &log_amax_arr, &Namax, 1);
	ReadArray("T_K.dat", &log_T_arr, &NT, 1);

	// compute spacing, necessary for interpolation
	// (assumes uniform spacing)
	d_q_arr = (q_arr[Nq-1] - q_arr[0]) / (Nq - 1);
	dlog_amax_arr = (log_amax_arr[Namax-1] - log_amax_arr[0]) / (Namax - 1);
	dlog_T_arr = (log_T_arr[NT-1] - log_T_arr[0]) / (NT - 1);

	// allocate memory for opacity arrays
	log_kR_arr = (double ***) Array3D(Nq, Namax, NT, sizeof(double));
	log_kP_arr = (double ***) Array3D(Nq, Namax, NT, sizeof(double));

	FILE *opac;
	int i,j,k;
	double tmp;

	// read opacity arrays from binary files, store as log values
	opac = fopen("kR_cm2g.dbl", "rb");

	for (k=0; k<Nq; k++) for (j=0; j<Namax; j++) for (i=0; i<NT; i++) {
		fread(&tmp, sizeof(double), 1, opac);
		log_kR_arr[k][j][i] = log(tmp);
	}

	fclose(opac);

	opac = fopen("kP_cm2g.dbl", "rb");

	for (k=0; k<Nq; k++) for (j=0; j<Namax; j++) for (i=0; i<NT; i++) {
		fread(&tmp, sizeof(double), 1, opac);
		log_kP_arr[k][j][i] = log(tmp);
	}

	fclose(opac);
}

double EvaluateRosselandOpacityArray(double q, double amax, double T) {
	/*
	Expects arguments in cgs: amax [cm], T [K].
	The function then uses amax in microns internally.
	This function computes the Rosseland mean opacity
	by trilinear interpolation in (q, log(amax), log(T)) space
	in the precomputed array.
	*/

	// convert amax to microns and take log for interpolation
	const double la = log(amax*1e4);
	const double lT = log(T);

	// find left indices for interpolation
	const int k = (int) ((Nq-1) * (q-q_arr[0]) / (q_arr[Nq-1] - q_arr[0]));
	const int j = (int) ((Namax-1) * (la-log_amax_arr[0]) / (log_amax_arr[Namax-1] - log_amax_arr[0]));
	const int i = (int) ((NT-1) * (lT-log_T_arr[0]) / (log_T_arr[NT-1] - log_T_arr[0]));

	// clamp indices to array bounds
	const int kL = MIN(MAX(k, 0), Nq-1);
	const int jL = MIN(MAX(j, 0), Namax-1);
	const int iL = MIN(MAX(i, 0), NT-1);

	// find right indices for interpolation. If at boundary, use same index -> no inter/extrapolation
	const int kR = kL == Nq-1    ? kL : kL+1;
	const int jR = jL == Namax-1 ? jL : jL+1;
	const int iR = iL == NT-1    ? iL : iL+1;

	// compute weights for interpolation
	const double wk = (q  - q_arr[kL]) / d_q_arr;
	const double wj = (la - log_amax_arr[jL]) / dlog_amax_arr;
	const double wi = (lT - log_T_arr[iL]) / dlog_T_arr;

	// indexing is q, a, T
	const double log_kappaLLL = log_kR_arr[kL][jL][iL], log_kappaLLR = log_kR_arr[kL][jL][iR];
	const double log_kappaLRL = log_kR_arr[kL][jR][iL], log_kappaLRR = log_kR_arr[kL][jR][iR];
	const double log_kappaRLL = log_kR_arr[kR][jL][iL], log_kappaRLR = log_kR_arr[kR][jL][iR];
	const double log_kappaRRL = log_kR_arr[kR][jR][iL], log_kappaRRR = log_kR_arr[kR][jR][iR];

	// collapse to index a, T using q weights
    const double log_kappaLL = log_kappaLLL * (1-wk) + log_kappaRLL * wk;
    const double log_kappaLR = log_kappaLLR * (1-wk) + log_kappaRLR * wk;
    const double log_kappaRL = log_kappaLRL * (1-wk) + log_kappaRRL * wk;
    const double log_kappaRR = log_kappaLRR * (1-wk) + log_kappaRRR * wk;

	// collapse to index T using amax weights
    const double log_kappaL  = log_kappaLL * (1-wj) + log_kappaRL * wj;
    const double log_kappaR  = log_kappaLR * (1-wj) + log_kappaRR * wj;

    // collapse to final value using T weights
    const double log_kappa = log_kappaL * (1-wi) + log_kappaR * wi;

	return exp(log_kappa);
}

double EvaluatePlanckOpacityArray(double q, double amax, double T) {
	/*
	Expects arguments in cgs: amax [cm], T [K].
	The function then uses amax in microns internally.
	This function computes the Planck mean opacity
	by trilinear interpolation in (q, log(amax), log(T)) space
	in the precomputed array.
	*/

	// convert amax to microns and take log for interpolation
	const double la = log(amax*1e4);
	const double lT = log(T);

	// find left indices for interpolation
	const int k = (int) ((Nq-1) * (q-q_arr[0]) / (q_arr[Nq-1] - q_arr[0]));
	const int j = (int) ((Namax-1) * (la-log_amax_arr[0]) / (log_amax_arr[Namax-1] - log_amax_arr[0]));
	const int i = (int) ((NT-1) * (lT-log_T_arr[0]) / (log_T_arr[NT-1] - log_T_arr[0]));

	// clamp indices to array bounds
	const int kL = MIN(MAX(k, 0), Nq-1);
	const int jL = MIN(MAX(j, 0), Namax-1);
	const int iL = MIN(MAX(i, 0), NT-1);

	// find right indices for interpolation. If at boundary, use same index -> no inter/extrapolation
	const int kR = kL == Nq-1    ? kL : kL+1;
	const int jR = jL == Namax-1 ? jL : jL+1;
	const int iR = iL == NT-1    ? iL : iL+1;

	// compute weights for interpolation
	const double wk = (q  - q_arr[kL]) / d_q_arr;
	const double wj = (la - log_amax_arr[jL]) / dlog_amax_arr;
	const double wi = (lT - log_T_arr[iL]) / dlog_T_arr;

	// indexing is q, a, T
	const double log_kappaLLL = log_kP_arr[kL][jL][iL], log_kappaLLR = log_kP_arr[kL][jL][iR];
	const double log_kappaLRL = log_kP_arr[kL][jR][iL], log_kappaLRR = log_kP_arr[kL][jR][iR];
	const double log_kappaRLL = log_kP_arr[kR][jL][iL], log_kappaRLR = log_kP_arr[kR][jL][iR];
	const double log_kappaRRL = log_kP_arr[kR][jR][iL], log_kappaRRR = log_kP_arr[kR][jR][iR];

	// collapse to index a, T using q weights
    const double log_kappaLL = log_kappaLLL * (1-wk) + log_kappaRLL * wk;
    const double log_kappaLR = log_kappaLLR * (1-wk) + log_kappaRLR * wk;
    const double log_kappaRL = log_kappaLRL * (1-wk) + log_kappaRRL * wk;
    const double log_kappaRR = log_kappaLRR * (1-wk) + log_kappaRRR * wk;

	// collapse to index T using amax weights
    const double log_kappaL  = log_kappaLL * (1-wj) + log_kappaRL * wj;
    const double log_kappaR  = log_kappaLR * (1-wj) + log_kappaRR * wj;

    // collapse to final value using T weights
    const double log_kappa = log_kappaL * (1-wi) + log_kappaR * wi;

	return exp(log_kappa);
}

int main() {
	// An example of how to use the opacity functions.
	const double q = -2.65;
	const double a = 0.4; // cm
	const double T = 120.0; // K
	// For this combination (-2.65, 0.4, 120) you should get:
	const double kR_exp = 5.191230176743109; // cm^2/g
	const double kP_exp = 4.102621899723951; // cm^2/g

	// You need to call this once
	ReadOpacityData();

	double kR, kP;
	const size_t start = clock();
	const int Ntries = 10000000;
	for (int i = 0; i < Ntries; i++) {
		// start timer
		// Now call opacity evaluation as many times as you'd want
		kR = EvaluateRosselandOpacityArray(q, a, T);
		kP = EvaluatePlanckOpacityArray(q, a, T);
		// stop timer
	}
	const size_t end = clock();
	const double duration = (double)(end - start) / CLOCKS_PER_SEC;
	
	printf("Duration: %.6f seconds\n", duration);
	printf("Average time per call: %.6f ns\n", duration / Ntries * 1e9);
	printf("kR = %.16le cm^2/g; error = %.16le\n", kR, kR - kR_exp);
	printf("kP = %.16le cm^2/g; error = %.16le\n", kP, kP - kP_exp);
	fflush(stdout);

	return 0;
}

// Make sure you have the data files in the same directory as this file:
// q.dat, amax_um.dat, T_K.dat, kR_cm2g.dbl, kP_cm2g.dbl

// compile this file (growpacity.c) into an executable called test with:
// gcc -o test growpacity.c -lm

// run with:
// ./test

// if it fails, you might need to adjust permissions:
// chmod +x test

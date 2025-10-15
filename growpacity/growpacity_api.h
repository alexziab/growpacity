#ifndef GROWPACITY_API_H
#define GROWPACITY_API_H

void ReadOpacityData();
double EvaluateRosselandOpacityArray(double q, double amax, double T);
double EvaluatePlanckOpacityArray(double q, double amax, double T);

#endif
#include <Python.h>
#include "growpacity_api.h"

static PyObject *py_EvaluateRosselandOpacityArray(PyObject *self, PyObject *args)
{
    double q, amax, T;
    if (!PyArg_ParseTuple(args, "ddd", &q, &amax, &T))
        return NULL;
    double result = EvaluateRosselandOpacityArray(q, amax, T);
    return PyFloat_FromDouble(result);
}

static PyObject *py_EvaluatePlanckOpacityArray(PyObject *self, PyObject *args)
{
    double q, amax, T;
    if (!PyArg_ParseTuple(args, "ddd", &q, &amax, &T))
        return NULL;
    double result = EvaluatePlanckOpacityArray(q, amax, T);
    return PyFloat_FromDouble(result);
}

static PyObject *py_ReadOpacityData(PyObject *self, PyObject *args)
{
    ReadOpacityData();
    Py_RETURN_NONE;
}

static PyMethodDef GrowpacityMethods[] = {
    {"read_opacity_data", py_ReadOpacityData, METH_NOARGS, "read the opacity data from file."},
    {"rosseland_opacity", py_EvaluateRosselandOpacityArray, METH_VARARGS, "compute the Rosseland opacity."},
    {"planck_opacity", py_EvaluatePlanckOpacityArray, METH_VARARGS, "compute the Planck opacity."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef growpacitymodule = {
    PyModuleDef_HEAD_INIT,
    "growpacity_c",
    NULL,
    -1,
    GrowpacityMethods};

PyMODINIT_FUNC PyInit_growpacity_c(void)
{
    return PyModule_Create(&growpacitymodule);
}
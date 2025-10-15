To test for coverage, make sure you have installed `pytest` and `pytest-cov`, then execute:

```NUMBA_DISABLE_JIT=1 pytest --cov growpacity  --cov-report html```

The `NUMBA_DISABLE_JIT` flag is required to include coverage for the interpolating function.

To test for coverage, make sure you have installed `pytest-cov`, then execute:

```NUMBA_DISABLE_JIT=1 pytest --cov growpacity  --cov-report html```

The `NUMBA_DISABLE_JIT` flag is required to include coverage for the interpolating function.

Coverage reports will be generated in the `htmlcov` folder. Open `htmlcov/index.html` in your browser to view the report.

The module achieves 96% coverage with both numba and scipy installed. The missing 4% corresponds to error handling during imports when either numba or scipy are not installed, which cannot be tested in the same environment. Nevertheless, these cases have been manually verified by uninstalling the respective packages and running the tests again with:
- without `scipy`
- without `numba`
- with neither `scipy` nor `numba`

For practical purposes, we have excluded these from automated coverage testing.
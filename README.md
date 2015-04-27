Poisson Disk Points Generator
=============================

This project wraps the [C++ Poisson Disk Generator by Sergey Kosarevsky](https://github.com/corporateshark/poisson-disk-generator) in Python using Cython. For examples see ``test_poisson.py``.

# Install

    python setup.py build_ext -i

# Run unit tests

```
py.test
py.test --interactive  # unskip interactive tests (matplotlib plots are
shown)
py.test --nocapturelog -s  # shows the log output
py.test test_poisson.py::TestPoisson::test_4_plot --interactive  # run a single test
```


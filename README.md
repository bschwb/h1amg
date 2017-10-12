# H1 Algebraic Multigrid for NGSolve

## Linux Build Steps
First make sure you have [NGSolve](https://ngsolve.org) installed.
You will also need [cmake](https://cmake.org/) and a C++ compiler.

After that clone this repository and build with cmake and make.
```
git clone https://github.com/bschwb/h1amg.git
mkdir -p h1amg/build
cd h1amg/build
cmake ../ && make install
```

## Examples
To run the python examples be sure to follow the build steps above.
Then navigate into the `python` subdirectory and run `netgen laplace_square.py`
or any other example in that directory.

## Testing
Tests are disabled by default.
To enable the tests invoke `cmake` with the flag `-DBUILD_TESTING=ON`.
To run the test navigate to the build directory and run `make test` or `ctest`.
If you need more output, e.g. to see specific tests failing, use `ctest -V`.
To run individual tests use `ctest -R <regex>`.
E.g. `ctest -R h1` to only run h1 integration tests.

## Profiling
At the beginning of the python script you want to profile add the following line to specify
the maximum tracefile size in bytes.
```python
import ngsolve as ngs
ngs.ngsglobals.pajetrace = 100000000
```

After running your example, a file named `ng0.trace` will be saved in the directory
you ran the example in.
You can view this file with [vite](http://vite.gforge.inria.fr/).

### Profiling Errors

If you get the following error message, try to increase the maximum tracefile size.
```
Tracing stopped during computation due to tracefile size limit of 0 megabytes.
To increase the limit, set in the pde file:
flags tracer = -max_size=size_in_megabytes

max_size=0 disables tracing
```

The error message describes how to fix it in case of using a `.pde` file which is not what we use
for here.

## Contributers

* [Joachim Schöberl](http://www.asc.tuwien.ac.at/~schoeberl/wiki/index.php/Joachim_Sch%C3%B6berl)
* Lukas Kogler
* Bernd Schwarzenbacher

## Linux Build Steps
This assumes installed Netgen/NGSolve.

```
git clone https://github.com/bschwb/h1amg.git
cd h1amg
mkdir build
cd build
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

## Contributers
Joachim Schöberl
Lukas Kogler
Bernd Schwarzenbacher

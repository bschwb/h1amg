add_python_test_target(h1_integration h1amg test_h1.py)
add_test(NAME h1_tests COMMAND make h1_integration)

Testing
=======

pytest
------

Basic testing is done with `pytest <https://docs.pytest.org/en/latest>`_.

Make sure that the code is installed from source as editable
```bash
python -m pip install -e ."[test]"
```

To test the code, go to where the source code is, and run
```bash
python -m pytest
```

tox
---

To test in different environments, you can use `tox <https://tox.readthedocs.io/en/latest/>`_. Install tox and run
```bash
tox
```

This tests all the environments listed in the `tox.ini` file, in `env_list`. To test a specific one, run e.g.
```bash
tox -e py39
```

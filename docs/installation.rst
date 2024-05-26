Installation
============

`popurri` is implemented in python 3.

Recommended: Create and activate an environment with e.g. conda and install `pip`:
```bash
conda create -n popurri
conda activate popurri
conda install pip
```

From source
-----------

The source code for `popurri` can be downloaded from GitHub and installed by running
```bash
git clone https://github.com/mlafarga/popurri.git
cd popurri
python -m pip install .
```

To install in editable mode and/or include also the dependencies for the documentation and testing, run:
```bash
python -m pip install -e .

# With docs dependencies
python -m pip install -e .[docs]

# With testing dependencies
python -m pip install -e .[test]

# With docs and testing dependencies
python -m pip install -e .[docs,test]
```

macOS note: zsh uses square brackets for globbing / pattern matching, so need to use e.g.
```bash
python -m pip install -e ."[docs,test]"
```

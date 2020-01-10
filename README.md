# SEIS glitch package (SEISglitch)

## Installation

Create a new ``conda`` environment and install dependencies

```
conda create -n ppol python=3.7 obspy pandas pyyaml -c conda-forge
conda activate ppol
```

Clone the repository

```
git clone https://gitlab.com/johnrobertscholz/ppol.git
cd ppol
```

Intall using ``pip``

```
pip install .
```

**or**

```
pip install -e .
```

for symbolic links to the module files, which avoids having to re-compile after editing the code. This is useful during code development.
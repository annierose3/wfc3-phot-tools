# wfc3-phot-tools

This Python library contains tools originally developed for the photometric calibration of the HST/WFC3 instrument.

This package is currently implemented in the IR staring mode, UVIS staring mode, and UVIS scanning mode monitoring pipelines. These monitors enable the WFC3 Instrument Team to evaluate the evolving sensitivity of the WFC3/IR and WFC3/UVIS detectors. As of June 2023, this library is being made publicly available in an effort to increase accessibility and transparency among our user base.

## Installation

### 1. Cloning the `wfc3-phot-tools` repository.

First, clone this repository and change into the newly created directory.

```
git clone https://github.com/spacetelescope/wfc3-phot-tools.git

cd wfc3-phot-tools
```

### 2. Creating a new environment (optional).

If desired, a new `conda` environment can be created with the necessary dependencies using the YAML file included in this repository, within which this package can be installed.

To do so, run the following command on the command line to create an environment called `wfc3_phot_env`.

```
conda env create -f wfc3_phot_env.yml
```

Then activate the virtual environment.

```
conda activate wfc3_phot_env
```

Alternately, this package can be installed as a standalone package, a process that will install the minimum required versions of the package's dependencies (listed in `pyproject.toml`) into the active environment.


### 3. Installing the `wfc3_phot_tools` library.

There are two ways to install this library. For most people, a static installation should suffice, as this is a stable package.

```
pip install .
```

Active developers of this package should install in "editable" mode.

```
pip install --editable .
```

## Contributors

The following people have contributed to this code base:
- Varun Bajaj
- Mariarosa Marinelli
- Jenny Medina
- Clare Shanahan

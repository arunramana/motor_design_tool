# dolomites package

![dolomites logo](logo.png)
## dolomites, a package for the design of electric machines and drives.

[https://gitlab.com/LuigiAlberti/dolomites-python](https://gitlab.com/LuigiAlberti/dolomites-python)

## Functionalities:
The various modules in the package perform the following computations
* koil:
  - design of balanced symmetrical windings
  - synthesis of custom windings
* apollo:
  from flux-linkages maps of synchronous machines to:
  - apparent and incremental inductances
  - mtpa trajectory
  - self-sensing capability evaluation
 * tiziano:
  some utilities to build FE meshes and run FE simulations
 * fnc:
  an utility to visualize space vectors related to 3-phase systems.
  Based on PySide6 widgets


## Installation on Windows
- Download and install the latest version of Python
- During the installation make sure to enable the option "Add Python to PATH"
- Download the source code "dolomites-python-master.zip"
- Extract the zip file
- Open the unzipped file
- Open the folder dolomites-python-master
- Launch the command prompt here and execute the following command:
```console
pip install .
```
- dolomites is now installed.
- Execute the following command to install Jupyter Notebook:
```console
pip install jupyter
```


## Installation on Linux
- python3 should be pre-installed on recent Linux distros.
- To install pip3 execute:
```console
sudo apt install python3-pip
```
- Download the source code "dolomites-python-master.zip"
- Extract the zip file
- Open the unzipped file
- Open the folder dolomites-python-master
- Open the terminal here and execute the following commands:
```console
pip3 install .
```
- dolomites is now installed.
- Execute the following command to install Jupyter Notebook:
```console
sudo apt install jupyter-notebook
```  

## Examples
Several examples are available in the tests directory.
Each module of dolomites has its own test directory: please explore them all!

Both Jupyter notebooks (.ipynb) and Python scripts (.py) are available.
For example, to run a Jupyter notebook launch a command similar to:
```console
jupyter-notebook apollo-test-linear-motor.ipynb
```  


## License and contributors
dolomites is copyright (C) 2006-2023 by L. Alberti and other contributors, University of Padova, and is distributed under the terms of the GNU General Public License (GPL) version 3.

Contributors are:

* Matteo Berto (apollo)
* Elia Scolaro (tiziano)
* Alice Maimeri (apollo)


# TALYS
TALYS is a software package for the simulation of nuclear reactions below 200 MeV. 
TALYS is based on state-of-art nuclear structure and reaction models. 

## Documentation and reference
A description of the code and its options can be found in the [TALYS Tutorial (pdf)](https://github.com/arjankoning1/talys/blob/main/doc/talys.pdf).
The reference to be used for TALYS is

Arjan Koning, Stephane Hilaire and Stephane Goriely, *TALYS: modeling of nuclear reactions*, European Journal of Physics A59 (6), 131 (2023).

## Installation

### Prerequisites:

The following are the prerequisites for compiling TALYS:
  - git (only if the package is downloaded via Github)
  - a recent Fortran compiler such as gcc (gfortran)

### Downloads:

To download TALYS, you can use one of the following options:
#### 1. Download the entire tar file (frozen version):
```
https://nds.iaea.org/talys/talys.tar
tar zxf talys.tar
```
#### 2. Using git (latest beta version):
```
git clone https://github.com/arjankoning1/talys.git
```
The TALYS structure database and sample cases do not fall under the git repository. Hence, to get a  working system you need to download
```
https://nds.iaea.org/talys/misc/structure.tar
https://nds.iaea.org/talys/samples/talys_samples.tar
```
and after
```
tar zxf structure.tar
tar zxf talys_samples.tar
```
you should move both *structure/* and *samples/* inside the *talys/* directory.

### Installation instructions :

To install TALYS, you can use one of the following options:
#### 1. Using make:
```
cd talys/source
make
```
#### 2. Using the install_talys.bash script:
```
cd talys
install_talys.bash talys
```

The above will produce a *talys* executable in the *talys/bin* directory. 
The compiler and its flags can be set in either *source/Makefile* or in *code_build.bash*.

### Memory restrictions:

For computers with (very) small RAM, or for installation on Windows, the total allocated memory may be too large. In that case, edit *A0_talys_mod.f90* in the source directory and reduce the value of the *memorypar* variable.

## The TALYS package

The *talys/* directory contains the following directories and files:

+ `README.md` this README file
+ `LICENSE` the License file
+ `install_talys.bash` and `code_build.bash` installation scripts
+ `source/` the Fortran source code of TALYS and the Makefile
+ `bin/` the *talys* executable after successful installation
+ `structure/` the nuclear structure and reaction database in various subdirectories
+ `misc/` miscellaneous files such as a gnuplot script to plot TALYS results versus EXFOR data
+ `doc/` the tutorial in pdf format
+ `samples/` the input and output files of the sample cases, and the *verify* script to run the sample cases

In total, you will need about 8 Gb of free disk space to install TALYS.

## Sample cases

The sample cases serve to provide examples of the use of TALYS and to verify a successful installation. The *samples/* directory contains various sample cases with a subdirectory *org/* with our results and a subdirectory *new/* with the results produced by the user. The entire sample set will take about 1 hour.
```
cd samples
./verify
```

You may create your own input file, e.g. *talys.inp* after which TALYS works as follows:
```
talys < talys.inp > talys.out
```
assuming that *talys* is linked to the *talys/bin/talys* executable.

## Plotting

The cross sections that TALYS calculates can be compared with experimental data from the EXFOR database using the misc/tplot script.
For this, the [EXFORTABLES](https://github.com/arjankoning1/exfortables) database needs to be installed in your home directory. 
Type 'misc/tplot' to get all the options for this plotting command. 'tplot' needs to be called from your working directory, 
i.e. the TALYS output files need to be present.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.

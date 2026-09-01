
# TALYS
TALYS is a software package for the simulation of nuclear reactions below 200 MeV. 
TALYS is based on state-of-the-art nuclear structure and reaction models. 

## Documentation and reference
A description of the code and its options can be found in the [TALYS Tutorial (pdf)](https://github.com/arjankoning1/talys/blob/main/doc/talys.pdf).
The reference to be used for TALYS is

Arjan Koning, Stephane Hilaire and Stephane Goriely, *TALYS: modeling of nuclear reactions*, European Physical Journal A59 (6), 131 (2023).

## Installation

### Prerequisites:

The following are the prerequisites for compiling TALYS:
  - git (only if the package is downloaded via GitHub)
  - GNU make
  - a recent Fortran compiler, such as GNU Fortran (gfortran)

### Downloads:

To download TALYS, you can use one of the following options:
#### 1. Download the entire tar file (frozen version TALYS-2.2):
This is available at the the [TALYS page](https://nds.iaea.org/talys/), and can be downloaded by clicking on the link or
```
curl -LO https://nds.iaea.org/talys/codes/talys.tar
tar zxf talys.tar
```
#### 2. Using git (latest beta version):
```
git clone https://github.com/arjankoning1/talys.git
```
### Installation instructions :

To install TALYS,
#### 1. For the tar file (frozen version TALYS-2.2)
```
cd talys
./install_talys.bash 
```
after which you will be prompted for your name, which will appear in the output files.
An alternative option is
```
cd talys/source
make
```
The above will invoke the default compiler gfortran.
The compiler and its flags can be set in either *source/Makefile* or in *code_build.bash*.
#### 2. For the git (latest beta) version,
```
cd talys
./install_talys.bash 
```
which automatically executes the *Makefile* in *talys/source*. At the end, *install_talys.bash*
will print the recommended shell configuration.
Set TALYS_DIR to the TALYS installation directory. This variable is
required unless the fallback path in source/machine.f90 has been set
manually. For example:
```
  export TALYS_DIR="/Users/koning/talys"
```
If you want to run *talys* from anywhere, set
```
  export PATH="$TALYS_DIR/bin:$PATH"
```
To include your name in the output files, set:
```
  export TALYS_USER="Your Name"
```

For the git version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, 
the Makefile will use the following flags:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

Compiler and compilation options can be passed to the Makefile through *install_talys.bash*.
e.g. you may replace the installation command above by
```
# GNU Fortran
./install_talys.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"

# Intel Fortran
./install_talys.bash FC=ifx FFLAGS="-O3"
```

The above will produce a *talys* executable in the *talys/bin* directory. 
### Memory restrictions:

For old computers with (very) small RAM, or for installation on Windows, the total allocated memory may be too large. 
In that case, edit *A0_talys_mod.f90* in the source directory and reduce the value of the *memorypar* variable.

## The TALYS package

The *talys/* directory contains the following directories and files:

+ `README.md` this README file
+ `LICENSE` the License file
+ `install_talys.bash` installation script
+ `source/` the Fortran source code of TALYS and the Makefile
+ `bin/` the *talys* executable after successful installation
+ `structure/` the nuclear structure and reaction database in various subdirectories
+ `misc/` miscellaneous files such as a gnuplot script to plot TALYS results versus EXFOR data
+ `doc/` the tutorial in pdf format
+ `samples/` the input and output files of the sample cases, and the *verify* script to run the sample cases

In total, you will need about 8 GB of free disk space to install TALYS.

## Sample cases

The sample cases serve to provide examples of the use of TALYS and to verify a successful installation. The *samples/* directory contains various sample cases with a subdirectory *org/* with our results and a subdirectory *new/* with the results produced by the user. The entire sample set will take about 1 hour.
```
cd samples
./verify
```
For the git version, you may do
```
make -C source check
```

You may create your own input file, e.g. *talys.inp* after which TALYS works as follows:
```
talys < talys.inp > talys.out
```
assuming that *talys/bin* has been added to PATH.

## Plotting

The cross sections that TALYS calculates can be compared with experimental data from the EXFOR database using the misc/tplot script.
For this, the [EXFORTABLES](https://github.com/arjankoning1/exfortables) database needs to be installed in your home directory. 
Type 'misc/tplot' to get all the options for this plotting command. 'tplot' needs to be called from your working directory, 
i.e. the TALYS output files need to be present.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.

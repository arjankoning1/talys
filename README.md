# TALYS

TALYS is a software package for the simulation of nuclear reactions below 200 MeV.  
TALYS is based on state-of-the-art nuclear structure and reaction models.

## Documentation and reference

A description of the code and its options can be found in the [TALYS Tutorial (pdf)](https://github.com/arjankoning1/talys/blob/main/doc/talys.pdf).

The reference to be used for TALYS is

Arjan Koning, Stephane Hilaire and Stephane Goriely, *TALYS: modeling of nuclear reactions*, European Physical Journal A59 (6), 131 (2023).

## Installation

### Prerequisites

The following are the prerequisites for compiling TALYS:

- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)
- git, only when TALYS is downloaded using `git clone`

### Downloads

TALYS can be downloaded in one of the following ways.

#### 1. Frozen version TALYS-2.2 (December 2025)

The frozen TALYS-2.2 distribution is available from the [TALYS page](https://nds.iaea.org/talys/). It can be retrieved by clicking on the download link or with

```bash
curl -LO https://nds.iaea.org/talys/codes/talys.tar
tar zxf talys.tar
```

This version is fixed and will not change.

#### 2. Latest beta version without git

Users who do not have git can download a snapshot of the current `main` branch directly from GitHub:

```bash
curl -L \
  -o talys-main.tar.gz \
  https://github.com/arjankoning1/talys/archive/refs/heads/main.tar.gz

tar zxf talys-main.tar.gz
mv talys-main talys
```

This produces the same `talys/` directory structure as the git version, but without the git history.

The downloaded snapshot contains the latest version of the `main` branch at the time of download. To obtain a newer version later, download the snapshot again.

#### 3. Latest beta version using git

Users with git can clone the repository with

```bash
git clone https://github.com/arjankoning1/talys.git
```

The advantage of this method is that the local TALYS installation can subsequently be updated with

```bash
cd talys
git pull --ff-only
```

### Installation instructions

#### 1. Frozen version TALYS-2.2

For the frozen tar distribution:

```bash
cd talys
./install_talys.bash
```

after which you will be prompted for your name, which will appear in the output files.

An alternative option is

```bash
cd talys/source
make
```

The above will invoke the default compiler gfortran.

The compiler and its flags can be set in either `source/Makefile` or `code_build.bash`.

#### 2. Latest beta version

The installation procedure is identical whether the latest beta version was obtained as a GitHub tar snapshot or using `git clone`.

From the `talys/` directory, run

```bash
./install_talys.bash
```

which automatically executes the Makefile in `talys/source`.

An alternative is

```bash
cd talys/source
make
```

The executable is installed as

```text
talys/bin/talys
```

At the end, `install_talys.bash` prints the recommended shell configuration.

Set `TALYS_DIR` to the TALYS installation directory. This variable is required unless the fallback path in `source/machine.f90` has been set manually. For example:

```bash
export TALYS_DIR="/Users/koning/talys"
```

If you want to run `talys` from anywhere, set

```bash
export PATH="$TALYS_DIR/bin:$PATH"
```

To include your name in the output files, set

```bash
export TALYS_USER="Your Name"
```

These lines can be added to `~/.zshrc` or `~/.profile`.

For the latest beta version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

Compiler and compilation options can be passed to the Makefile through `install_talys.bash`. For example:

```bash
# GNU Fortran
./install_talys.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"

# Intel Fortran
./install_talys.bash FC=ifx FFLAGS="-O3"
```


### Memory restrictions

For old computers with very small RAM, or for installation on Windows, the total allocated memory may be too large.

In that case, edit `A0_talys_mod.f90` in the source directory and reduce the value of the `memorypar` variable.

## The TALYS package

The `talys/` directory contains the following directories and files:

- `README.md` this README file
- `LICENSE` the license file
- `install_talys.bash` installation script
- `source/` the Fortran source code of TALYS and the Makefile
- `bin/` the `talys` executable after successful installation
- `structure/` the nuclear structure and reaction database in various subdirectories
- `misc/` miscellaneous files such as a gnuplot script to plot TALYS results versus EXFOR data
- `doc/` the tutorial in pdf format
- `samples/` the input and output files of the sample cases and the `verify` script to run the sample cases

In total, you will need about 8 GB of free disk space to install TALYS.

## Sample cases

The sample cases provide examples of the use of TALYS and can be used to verify a successful installation.

The `samples/` directory contains various sample cases with a subdirectory `org/` containing the reference results and a subdirectory `new/` containing results produced by the user.

The entire sample set takes about 1 hour.

To run the sample cases:

```bash
cd samples
./verify
```

For the latest beta version, the same check can be started from the top-level TALYS directory with

```bash
make -C source check
```

You may create your own input file, for example `talys.inp`, after which TALYS can be run as

```bash
talys < talys.inp > talys.out
```

assuming that `talys/bin` has been added to `PATH`.

## Plotting

The cross sections calculated by TALYS can be compared with experimental data from the EXFOR database using the `misc/tplot` script.

For this, the [EXFORTABLES](https://github.com/arjankoning1/exfortables) database needs to be installed in your home directory.

Run

```bash
misc/tplot
```

to see all available options.

`tplot` needs to be called from the working directory containing the TALYS output files.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.

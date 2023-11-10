Copyright (C) 2002 Regents of the University of Michigan,
portions used with permission.
For more information, see http://csem.engin.umich.edu/tools/swmf

This document outlines how to install the SWMF on your system and how
to create and access the documentation.  To learn more about the SWMF,
including how to compile and run the code, please consult the user
manual.  To install the SWMF and create the user manual please follow
the instructions below.

# Obtain SWMF

Get the full source code from GitHub/SWMFsoftware or the open-source code from
GitHub/MSTEM-QUDA.

The minimum requirement is the `SWMF` repository. 

You may also need the open-source `SWMF_data` repository that contains
large data files and can be downloaded into your home directory (or into
the SWMF directory):

```
cd
git clone https://github.com/SWMFsoftware/SWMF_data --depth=1
```

For solar applications (solar corona, inner heliosphere, CMEs) the
`SWMFSOLAR` repository can be useful. This is typically downloaded
into the SWMF directory, or to the (scratch/nobackup...) disk of
a supercomputer where the runs are performed:
```
cd
git clone https://github.com/SWMFsoftware/SWMFSOLAR --depth=1
```

Some data files used by the Center for Radiative Shock Hydrodynamics (CRASH)
are in the `CRASH_data` repository that is available to registered users.
If needed, it has to be placed into the home directory. 

# Getting the open-source MSTEM-QUDA/SWMF:

Clone the SWMF from GitHub/MSTEM-QUDA

```
cd {where_you_want_to_have_mstem-quda}
git clone https://github.com/MSTEM-QUDA/SWMF
```

The rest of the repositories (share, util, BATSRUS ...)
will be cloned from GitHub/MSTEM-QUDA during the installation.

# Getting the full SWMF from GitHub/SWMFsoftware (requires access)

Read the
[Git instructions](http://herot.engin.umich.edu/~gtoth/SWMF/doc/Git_instructions.pdf)
about registering, passwordless access, mail notifications, and
using the [gitclone](https://github.com/SWMFsoftware/share/blob/master/Scripts/gitclone) script.

## Clone the SWMF repository

```
cd {where_you_want_to_have_the_swmf}
gitclone SWMF
```

## Clone the CRASH_data repository into the home directory if needed
```
cd
gitclone CRASH_data
```

# Install SWMF

Many machines used by UofM are already recognized by the
`share/Scripts/Config.pl` script, which is called by all other `Config.pl`
scripts in the SWMF.
For these platform/compiler combinations installation is very simple:
```
./Config.pl -install
```
On other platforms the Fortran (and C) compilers should be explicitly given.
To see available choices, type
```
./Config.pl -compiler
```
Then install the code with the selected Fortran (and default C) compiler, e.g.
```
./Config.pl -install -compiler=gfortran
```
A non-default C compiler can be added after a comma, e.g.
```
./Config.pl -install -compiler=mpxlf90,mpxlc
```
For machines with no MPI library, use
```
./Config.pl -install -nompi -compiler=....
```
This will only allow serial execution, of course. Like with most scripts
in the SWMF, type
```
./Config.pl -help
```
for a comprehensive description of options and examples.

The ifort compiler (and possibly others too) use the stack for temporary
arrays, so the stack size should be large. For csh/tcsh add the following
to `.cshrc`:
```
unlimit stacksize
```
For bash/ksh/zsh add the following to `.bashrc` or equivalent initialization
file:
```
ulimit -s unlimited
```

# Create the manuals

Please note that creating the PDF manuals requires that LaTeX
(available through the command line) is installed on your system.

To create the PDF manuals type
```
make PDF
```
in the SWMF directory. The manuals will be in the doc/ directories.

## Cleaning the documentation
```
cd doc/Tex
make clean
```
To remove all the created documentation type
```
cd doc/Tex
make cleanpdf
```
As for most Makefile-s in the SWMF, type
```
make help
```
for a comprehensive list of make targets and examples.

# Read the manuals

All manuals can be accessed by opening the top index file 
```
open doc/index.html
```
You may also read the PDF files directly with a PDF reader.
The most important document is the user manual in
```
doc/SWMF.pdf
```

# Running tests

You can try running the standard test suite by typing
```
make -j test
```
in the main directory. The `-j` flag allows parallel compilation.
This requires a machine where `mpiexec` is available.
The tests run on 2 CPU cores by default.
The results of the tests are summarized in `test_swmf.res`.
Successful passing of the test is indicated by empty `*.diff` files.

To run the tests on more (up to 8) cores use
```
make -j test NP=4
```
You can also run an individual test. The list of available SWMF tests can be listed with
```
make test_help
```
For example, to run test1 without MPI on a single core use
```
make -j test1 MPIRUN=
```

SWMF
====

This document outlines how to install the SWMF on your system and how to
create and access the documentation.  To learn more about the SWMF, including
how to compile and run the code, please consult the user manual.  To install
the SWMF and create the user manual please follow the instructions below. 

Obtain SWMF
-----------

Get the source code from Git or from the tar balls.

The minimim requirement is the SWMF repository. 

You may also need the `SWMF_data` repository that contains large data files.
It can be put into the home directory or into the SWMF directory.

Some data files used by the Center for Radiative Shock Hydrodynamics (CRASH)
are in the `CRASH_data` repository.
If needed, it has to be placed into the home directory.

Getting SWMF from Git (requires access to UM GitLab)
----------------------------------------------------

Read the [GitLab instructions](http://herot.engin.umich.edu/~gtoth/SWMF/doc/GitLab_instructions.pdf) about registering, passwordless access, mail notifications, and defining the "gitlabclone" (or gitclone) alias/function.  

Check out the SWMF distribution

```bash
cd {where_you_want_to_have_the_swmf}
gitlabclone SWMF
```

Check out the `SWMF_data` distribution into the home directory if needed

```bash
cd
gitlabclone SWMF_data
```

Check out the `CRASH_data` distribution into the home directory if needed

```bash
cd
gitlabclone CRASH_data
```

### Opening tar balls

Open SWMF tar ball

```bash
cd {where_you_want_to_have_the_swmf}
tar -xzf {path_to_file}/SWMF_v{version_number}.tgz
```

If needed, open the `SWMF_data` tar ball into the home directory

```bash
cd
tar -xzf {path_to_file}/SWMF_data.tgz
```

If needed, open the `CRASH_data` tar ball into the home directory

```bash
cd
tar -xzf {path_to_file}/CRASH_data.tgz
```

Install SWMF
------------

Many machines used by UofM are already recognized by the share/Scripts/Config.pl script which is called by all other Config.pl scripts in the SWMF. For these platform/compiler combinations installation is very simple:

```bash
Config.pl -install
```

On other platforms the Fortran (and C) compilers should be explicitly given.  To see available choices, type 

```bash
Config.pl -compiler
```

Then install the code with the selected Fortran (and default C) compiler, e.g.

```bash
Config.pl -install -compiler=gfortran
```

A non-default C compiler can be added after a comma, e.g.

```bash
Config.pl -install -compiler=mpxlf90,mpxlc
```

For machines with no MPI library, use

```bash
Config.pl -install -nompi -compiler=....
```

This will only allow serial execution, of course.

The ifort compiler (and possibly others too) use the stack for temporary arrays,
so the stack size should be large. For csh/tcsh add the following to .cshrc: 

```bash
unlimit stacksize
```

For bash/ksh add the following to .bashrc or equivalent initialization file:

```bash
ulimit -s unlimited
```

Create the manuals
------------------

Please note that creating the PDF manuals requires that LaTex (available through the command line) and ps2pdf be installed on your system.

To create the PDF manuals type

```bash
make PDF
```

in the SWMF directory. The manuals will be in the doc/ directory.

Cleaning the documentation

```bash
cd doc/Tex
make clean
```

To remove all the created documentation type

```bash
cd doc/Tex
make cleanpdf
```

Read the manuals
----------------

All manuals can be accessed by opening the top index file 

```bash
open doc/index.html
```

You may also read the PDF files directly with a PDF reader.  The most important document is the user manual in

```bash
doc/SWMF.pdf
```

Running tests
-------------

You can try running the standard test suite by typing

```bash
make -j test
```

in the main directory. The -j flag allows parallel compilation.  This requires a machine where mpirun is available.  The tests run on 2 CPU cores by default.  The results of the tests are summarized in `test_swmf.res` Successful passing of the test is indicated by empty *.diff* files.

To run the tests on more (up to 8) cores use

```bash
make -j test NP=4
```

You can also run an individual test. The list of available SWMF tests can be listed with

```bash
make test_help
```

For example, to run test3 without MPI on a single core use

```bash
make -j test3 MPIRUN=
```

Copyright (C) 2002 Regents of the University of Michigan, portions used with permission. For more information, see http://csem.engin.umich.edu/tools/swmf


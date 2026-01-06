This program aims to search for a water molecule's energy positions in cavities of proteins having minimum energy. For more information, please refer to the original [Zhang and Hermans](https://www.ncbi.nlm.nih.gov/pubmed/9162944) paper.

Dowser's official link is (was?) http://danger.med.unc.edu/hermans/dowser/dowser.htm. However, it is hardly ever available, which is why I decided to upload the Dowser code to ensure program availability.

## Introduction

Prior to installing Dowser, update the system package index. For example, on Debian based distributions, type:

    sudo apt update
    sudo apt install build-essential

You need a FORTRAN (e.g. `gfortran`, `ifort`) and C (e.g. `gcc`) compilers in Linux to install Dowser. These can be installed on Debian based distributions with:

    sudo apt-get install gfortran gcc

On this repository, the make configuration in file `CODE/Makearch.linux` already includes `gfortran` as the FORTRAN compiler. If you are using another compiler, or in case you are using any other of the available architectures, please modify the line `F77		= f77` inside the corresponding make configuration file to meet your specific situation.

Also, this code was written back in the day using the K&R C (pre-ANSI) standard. An attempt to compile it with a modern C standard results in compilation errors due to old-style function definitions. To ignore them, you need to include the flags `-Wno-old-style-definition -std=gnu89` into the `ARCH_COPTS` list of your architecture make file.

If you are using the Linux architecture, this has already been done for you (see `CODE/Makearch.linux`), otherwise, include them by yourself into the equivalent `CODE/Makearch.[your-architecture]` before attempting to compile the Dowser source code.

## Installation

Once you have downloaded the code, get into its directory, and grant the following scripts with execution rights if they don't already have it. Then execute the installation:

     chmod u+x Install bin/dowser bin/append_dowser_water bin/append_index bin/cleanup_dowser bin/dryer bin/ms2pdb
     ./Install

As a result, you should obtain a series of new executable under directory `bin/linux/` (or any other architecture you had selected). Make sure to grant them execution rights too:

    chmod u+x bin/linux/*

Now, to make these executable accessible from anywhere on the system, export the Dowser program location to your PATH. For BASH, edit your .bashrc or .bash_profile file by adding:

     export DOWSER=/path/to/dowser  # write your Dowser absolute installation path here
     export DOW_MACH=linux
     export PATH=$PATH:$DOWSER/bin:$DOWSER/bin/$DOW_MACH

For CSH type:

     setenv DOWSER /path/to/dowser  # write the Dowser absolute installation path here
     setenv DOW_MACH linux
     set path = ( $path $DOWSER/bin $DOWSER/bin/$DOW_MACH )

And finally, source the edited file (whichever is your case) to make the changes take effect:

     source .bashrc

or

     source .bash_profile
    

## Residues in `DATA/atomdict.db`

To meet [pyARMm](https://github.com/yaideltg/pyarmm) requirements, new residues have been added to file `DATA/atomdict.db`. Furthermore, it needs to be pointed out that residue `HIS TERM NH3 COO` has been modified from

```
RESIDUE HIS TERM NH3 COO
ATOM HIS  N    C    CA    1.320  114.0  180.0  -0.280 N
ATOM HIS  H    N    NOT   1.000  123.0    0.0   0.280 H
ATOM HIS  CA   N    C     1.470  123.0  180.0   0.000 CH1
ATOM HIS  CB   CA   CG    1.530  110.0   60.0   0.000 CH2
ATOM HIS  CG   CB   ND1   1.530  112.0  180.0   0.000 CR
ATOM HIS  CD2  CG   NOT   1.360  132.0  180.0   0.000 CHR
ATOM HIS  ND1  CG   CE1   1.390  122.0    0.0   0.128 N
ATOM HIS  HD1  ND1  NOT   1.000  126.0     0.0  0.192 H
ATOM HIS  CE1  ND1  NE2   1.320  108.0  180.0   0.259 CHR
ATOM HIS  NE2  CE1  NOT   1.310  109.0    0.0  -0.579 N
ATOM HIS  C    CA   N     1.530  110.0  180.0   0.380 CR
ATOM HIS  O    C    NOT   1.240  121.0    0.0  -0.380 O
```

to

```
RESIDUE HIS TERM NH3 COO
ATOM HIS  N    C    CA    1.320  114.0  180.0  -0.3479 N
ATOM HIS  H    N    NOT   1.000  123.0    0.0   0.2747 H
ATOM HIS  CA   N    C     1.470  123.0  180.0  -0.1354 CH1
ATOM HIS  CB   CA   CG    1.530  110.0   60.0  -0.0414 CH2
ATOM HIS  CG   CB   ND1   1.530  112.0  180.0  -0.0012 CR
ATOM HIS  CD2  CG   NOT   1.360  132.0  180.0  -0.1141 CHR
ATOM HIS  ND1  CG   CE1   1.390  122.0    0.0  -0.1531 N
ATOM HIS  HD1  ND1  NOT   1.000  126.0    0.0   0.2317 H
ATOM HIS  CE1  ND1  NE2   1.320  108.0  180.0  -0.0170 CHR
ATOM HIS  NE2  CE1  NOT   1.310  109.0    0.0  -0.1718 N
ATOM HIS  HE2  NE2  NOT   1.000  126.0  180.0   0.2681 H
ATOM HIS  C    CA   N     1.530  110.0  180.0   0.7341 CR
ATOM HIS  O    C    NOT   1.240  121.0    0.0  -0.5894 O
```


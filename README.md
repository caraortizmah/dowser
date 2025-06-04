This program aims to search for a water molecule's energy positions in cavities of proteins having minimum energy. For more information, please refer to the original [Zhang and Hermans](https://www.ncbi.nlm.nih.gov/pubmed/9162944) paper.

Dowser's official link is (was?) http://danger.med.unc.edu/hermans/dowser/dowser.htm. However, it is hardly ever available, which is why I decided to upload the Dowser code to ensure program availability.

## Introduction

Prior to installing Dowser, update the system package index. For example, on Debian based distributions, type:

    sudo apt update
    sudo apt install build-essential

You need a FORTRAN (e.g. `gfortran`, `ifort`) and C (e.g. `gcc`) compilers in Linux to install Dowser. These can be installed on Debian based distributions with:

    sudo apt-get install gfortran gcc

On this repository, the make configuration in file `CODE/Makearch.linux` already includes `gfortran` as the FORTRAN compiler. If you are using another compiler, or in case you are using any other of the available architectures, please modify the line `F77		= f77` inside the corresponding make configuration file to meet your specific situation.

## Installation

Once you have downloaded the code get into its folder, and execute the following orders:

     chmod u+x Install
     chmod u+x bin/dowser
     ./Install
     
When installation ends, you need to export the Dowser program location to your PATH. For BASH, edit your .bashrc or .bash_profile file by typing:

     export DOWSER=   # write the Dowser absolute installation path here
     export DOW_MACH=linux
     export PATH=$PATH:$DOWSER/bin:$DOWSER/bin/$DOW_MACH

For CSH type:

     setenv DOWSER    # write the Dowser absolute installation path here
     setenv DOW_MACH linux
     set path = ( $path $DOWSER/bin $DOWSER/bin/$DOW_MACH )

And then, source the edited file (whichever in your case):

     source .bashrc

or

     source .bash_profile
    
That's all there is to it.

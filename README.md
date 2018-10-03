# What is polyDFE?
polyDFE infers the distribution of fitness effects (DFE) of new mutations
from polymorphism, and, if available, divergence data obtained from an
outgroup. The polymorphism data is given through the unfolded site frequency
spectrum (SFS). Note that, in order to obtain an unfolded SFS, at least one
outgroup is needed. However, including the divergence data when inferring
the DFE can lead to bias in the estimated parameters.

Once the DFE is obtained, polyDFE also calculates alpha, the rate of adaptive
molecular evolution.

polyDFEv2.0 can fit jointly multiple datasets for which some parameters can
be set to be invariant (shared across the datasets). This enables model testing
for parameter invariance, for example, testing if different datasets have the
same DFE or not.

More details can be found in
1. Tataru P., Mollion M., Glemin S., Bataillon T. (2017). Inference of
distribution of fitness effects and proportion of adaptive substitutions from
polymorphism data. Genetics, 207(3):1103-1119
1. Tataru P., Bataillon T. (2018). polyDFEv2.0: Testing for invariance of the
distribution of fitness effects within and across species. bioRxiv,
doi: https://doi.org/10.1101/363887

If using this software in a publication, please cite the above articles.


# Contact
Paula Tataru, paula.tataru at gmail . com


# Requirements
polyDFE can be ran using the provided pre-compiled binaries for Windows, Linux and iOS.

polyDFE is implemented in C and uses the GSL library. For compiling polyDFE,
the GSL library should be downloaded from https://www.gnu.org/software/gsl/ and
compiled accordingly. polyDFE has been tested using GSL-1.16. Note that polyDFE
does not work with a newer version of GSL.


# How to install polyDFE
The simplest way to compile polyDFE:

1. `cd` to the directory containing the makefile and type
`make all` to compile the code.

1. You can remove the program binaries and object files from the
source code directory by typing `make clean`.

If the GSL library is not installed in the default directory, than the
installation path needs to be specified in the makefile by updating USR_INC and
USR_LIB and uncommenting the lines at the top of the makefile.

polyDFEv2.0 is also distributed as pre-compiled binaries for Windows, Linux and macOS.


# How to run polyDFE
For more details on running polyDFE, please consult the manual and tutorial.

# License information
See the file LICENSE.txt for information on terms & conditions for usage,
and a DISCLAIMER OF ALL WARRANTIES.

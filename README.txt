What is polyDFE?
============================
polyDFE infers the distribution of fitness effects (DFE) of new mutations 
from polymorphism, and, if available, divergence data obtained from an 
outgroup. The polymorphism data is given through the unfolded site frequency 
spectrum (SFS). Note that, in order to obtain an unfolded SFS, at least one 
outgroup is needed. However, including the divergence data when inferring 
the DFE can lead to bias in the estimated parameters. 

Once the DFE is obtained, polyDFE also calculates alpha, the rate of adaptive 
molecular evolution.

More details can be found in 
[1] Tataru P., Mollion M., Glemin S., Bataillon T. (2016). Inference of 
distribution of fitness effects and proportion of adaptive substitutions from 
polymorphism data. Genetics, submitted.
[2] Tataru P., Mollion M., Franco M.A.P., Bataillon T. (2016). polyDFE: 
inference of distribution of fitness effects and proportion of adaptive 
substitutions. Bioinformatics, in preparation.

If using this software in a publication, please cite the above articles.


Contact
============================
Paula Tataru, paula@birc.au.dk


Requirements
============================
polyDFE is implemented in C and uses the GSL library. This should be downloaded
from https://www.gnu.org/software/gsl/ and compiled accordingly. polyDFE has 
been tested using GSL-1.16.


How to install polyDFE
============================
The simplest way to compile polyDFE:

  1. `cd' to the directory containing the makefile and type
     `make all' to compile the code.

  2. You can remove the program binaries and object files from the
     source code directory by typing `make clean'.
     
If the GSL library is not installed in the default directory, than the
installation path needs to be specified in the makefile by updating USR_INC and 
USR_LIB and uncommenting the lines at the top of the makefile.


How to run polyDFE
============================
For more details on running polyDFE, please consult the manual.

License information
============================
See the file LICENSE.txt for information on terms & conditions for usage,
and a DISCLAIMER OF ALL WARRANTIES.

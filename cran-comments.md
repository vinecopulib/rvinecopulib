## Breaking changes

This update breaks the current version of the 'vinereg' package. This is caused
by changes in the internal representation of models in our library. 'vinereg' 
modifies corresponding S3 objects manually in a way that 
is no longer valid. I am also the maintainer of the vinereg package and will
submit a patched version as soon our update is on it's way to CRAN. 

## Non-standard version number

The version (0.3.0.1.0) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* macOS High Sierra 10.13.6 (local install), R 3.5.0 with gcc-8
* ubuntu 14.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), R 3.5.0

## R CMD check results
There were no ERRORs or WARNINGs. 

## Reverse dependencies
Checked simIReff: 0 errors | 0 warnings | 0 notes
Checked vinereg : 2 errors | 2 warnings | 1 note

## Breaking changes

This update breaks the current version of the 'vinereg' package. This is caused
by changes in the internal representation of models in our library. 'vinereg' 
modifies corresponding S3 objects of manually in a way that 
is no longer valid. I am also the maintainer of the vinereg package and will
submit a patched version as soon our update is accepted. 

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
There were no ERROR, WARNINGs or NOTES. 

## Reverse dependencies

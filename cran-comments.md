## Non-standard version number
The version (0.2.4.1.0) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* local OS X install, R 3.4.2, using both gcc and apple's LLVM
* ubuntu 12.04 (on travis-ci), R 3.4.1
* Windows Server 2012 R2 x64 and x86 (on appveyor), R 3.4.1

## R CMD check results
There were no ERROR or WARNINGs. 

## Reverse dependencies
There are no reverse dependencies.

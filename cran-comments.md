## Non-standard version number
The version (0.2.5.1.0) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* macOS Sierra 10.12.6 (local install), R 3.4.2 with both Apple clang and GNU gcc
* ubuntu 14.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), R 3.4.1

## R CMD check results
There were no ERROR, WARNINGs or NOTES. 

## Reverse dependencies
There are no reverse dependencies.

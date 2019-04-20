## Non-standard version number

The version (0.3.1.1.0) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* ubuntu 14.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), R 3.5.3

## R CMD check results
There were no ERRORs or WARNINGs. 
There is a NOTE related to the size of the library, following the removal of the unconditional stripping of the library.

## Reverse dependencies
Checked simIReff: 0 errors | 0 warnings | 0 notes
Checked vinereg : 0 errors | 0 warnings | 0 notes

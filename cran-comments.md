## Non-standard version number

The version (0.5.0.1.0) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* ubuntu 16.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), devel

## R CMD check results
There were no ERRORs or WARNINGs. 

## Reverse dependencies
Checked simIReff: 0 errors | 0 warnings | 0 notes
Checked vinereg : 0 errors | 0 warnings | 0 notes

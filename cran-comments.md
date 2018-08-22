Patch following an email by Prof. Ripley warning us that the Solaris build
failed due to non-portable use of `log()` in C++. The patched version uses
`std::log()`, which is portable.

## Non-standard version number

The version (0.3.0.1.1) does not have three components. Since the package is 
an interface to the C++ library vinecopulib, we decided to adopt the convention 
of e.g. RcppEigen, where the first three components are related to 
vinecopulib's version.

## Test environments
* macOS High Sierra 10.13.6 (local install), R 3.5.0 with gcc-8
* ubuntu 14.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), R 3.5.1

## R CMD check results
There were no ERRORs or WARNINGs. 

## Reverse dependencies
Checked simIReff: 0 errors | 0 warnings | 0 notes
Checked vinereg : 0 errors | 0 warnings | 0 notes

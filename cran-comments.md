Urgent patch fixing a heap-buffer-overflow detected by the address sanitizer:
https://cran.r-project.org/web/checks/check_results_rvinecopulib.html

## Test environments
* ubuntu 16.04 (travis-ci), release/oldrel/devel
* Windows Server 2012 R2 x64 and x86 (appveyor), devel
* ubuntu 18.04 with gcc ASAN (release)
* ubuntu 18.04 with valgrind (release)

## R CMD check results
There were no ERRORs or WARNINGs. 

## Reverse dependencies
Checked simIReff: 0 errors | 0 warnings | 0 notes
Checked vinereg : 0 errors | 0 warnings | 0 notes

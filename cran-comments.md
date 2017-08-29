## Test environments
* local OS X install, R 3.4.1
* ubuntu 12.04 (on travis-ci), R 3.4.1
* Windows Server 2012 R2 x64 and x86 (on appveyor), R 3.4.1

## R CMD check results
There were no ERROR or WARNINGs. 
There was a NOTE related to sub-directories of 1Mb or more (only for the 
ubuntu builts). It was reproduced on a local ubuntu 12.04 install but 
everything appears to be in order.

## Reverse dependencies
There are no reverse dependencies.

Update following an email by Brian Ripley. The library is now compiled without 
forcing c++11, but boost uses features that were deprecated in later versions. 
As suggested by Ripley, we now add the `-D_HAS_AUTO_PTR_ETC=0` include flag to
disable the problematic parts.

## Test environments

* CRAN win builder (release)
* Windows Server 2019 (release)
* macOS (release)
* ubuntu 20.04 (release, oldrel, devel)


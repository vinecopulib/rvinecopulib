#!/bin/bash

git clone --depth 1 https://github.com/vinecopulib/vinecopulib/ -b prepare-v0.3.0 --single-branch

rm -rf ../inst/include/vinecopulib*
mv ./vinecopulib/include/* ../inst/include

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

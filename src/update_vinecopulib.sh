#!/bin/bash

git clone --depth 1 git@github.com:vinecopulib/vinecopulib.git -b dev --single-branch

rm -rf ../inst/include/vinecopulib*
mv ./vinecopulib/include/* ../inst/include

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

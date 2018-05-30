#!/bin/bash

git clone https://github.com/vinecopulib/vinecopulib/
cd vinecopulib
git checkout thread-includes
cd ..

rm -rf ../inst/include/vinecopulib*
mv ./vinecopulib/include/* ../inst/include

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

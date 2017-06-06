#!/bin/bash

git clone git@github.com:vinecopulib/vinecopulib.git

cd vinecopulib 
git checkout dev
cd ..

rm -rf bicop 
rm -rf misc 
rm -rf vinecop
rm -rf ../inst/include/vinecopulib*

mv ./vinecopulib/src/* .
mv ./vinecopulib/include/* ../inst/include

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp
sed -i '12i#include <RcppEigen.h>' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

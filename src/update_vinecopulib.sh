#!/bin/bash

git clone git@github.com:vinecopulib/vinecopulib.git

rm -rf bicop 
rm -rf misc 
rm -rf vinecop
rm -rf ../inst/include

mv ./vinecopulib/src/* .
mv ./vinecopulib/include ../inst

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp
sed -i '12i#include <RcppEigen.h>' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

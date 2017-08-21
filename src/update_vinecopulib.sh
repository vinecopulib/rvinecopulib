#!/bin/bash

git clone https://github.com/vinecopulib/vinecopulib
cd vinecopulib
git checkout tlls
cd ..

rm -rf bicop
rm -rf misc
rm -rf vinecop
rm -rf ../inst/include/vinecopulib*

mv ./vinecopulib/src/* .
mv ./vinecopulib/include/* ../inst/include

gsed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp
gsed -i '12i#include <RcppEigen.h>' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/mainpage.h
rm -rf vinecopulib

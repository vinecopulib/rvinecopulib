#!/bin/bash

git clone --depth 1 git@github.com:vinecopulib/vinecopulib.git -b discrete-fix  --single-branch

rm -rf ../inst/include/vinecopulib/*
mv ./vinecopulib/include/* ../inst/include

sed -i '11i#define INTERFACED_FROM_R' ./../inst/include/vinecopulib/misc/tools_interface.hpp

rm -rf ./../inst/include/vinecopulib/mainpage.h
rm -rf vinecopulib

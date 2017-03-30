#!/bin/bash

git clone git@github.com:vinecopulib/vinecopulib.git

rm -rf bicop 
rm -rf misc 
rm -rf vinecop
rm -rf ../inst/include
mv ./vinecopulib/src/* .
mv ./vinecopulib/include ../inst

rm -rf vinecopulib

#!/usr/bin/env bash

CWD=$(pwd)
cd ../doc
doxygen doxyfile
cd $CWD

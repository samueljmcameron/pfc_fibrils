#!/bin/bash

make

./EvsR $1 $2 $3

(cd ../Results; python plot.py $1 $2 $3 $4)
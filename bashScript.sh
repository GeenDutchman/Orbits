#!/bin/bash
echo starting script


python3 binaryNewton.py > binary1.txt

gnuplot 'binPlot.txt'



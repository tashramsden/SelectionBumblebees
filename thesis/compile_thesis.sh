#!/bin/bash
# Author: Tash Ramsden ter21@imperial.ac.uk
# Script: compile_thesis.sh
# Description: Bash script to compile LaTeX

pdflatex BumblebeeSelection.tex
bibtex BumblebeeSelection
pdflatex BumblebeeSelection.tex
pdflatex BumblebeeSelection.tex

rm *.aux
rm *.log
rm *.bbl
rm *.blg
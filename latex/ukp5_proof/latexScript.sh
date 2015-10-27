#!/bin/bash

pdflatex "$1"
bibtex "${1%.*}"
pdflatex "$1"
pdflatex "$1"
rm "${1%.*}.aux"


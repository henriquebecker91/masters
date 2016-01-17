#!/bin/bash

tex='ukp5.tex'
pdflatex $tex
bibtex "${tex%.*}.bib"
pdflatex $tex
pdflatex $tex
rm "${tex%.*}.aux"


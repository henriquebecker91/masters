#!/bin/bash

tex='slides_sea2016.tex'
pdflatex "$tex"
# Execute a second time to create the outline
pdflatex "$tex"
rm *.aux *.log *.toc *.nav *.out *.snm


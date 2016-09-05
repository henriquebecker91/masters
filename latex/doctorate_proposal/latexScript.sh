#!/bin/bash

texfile='plano-de-trabalho.tex'
pdflatex "$texfile"
bibtex "${texfile%.*}"
pdflatex "$texfile"
pdflatex "$texfile"
#rm "${texfile%.*}.aux"


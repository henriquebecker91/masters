#!/bin/bash

tex='sbc-template.tex'
pdflatex $tex
bibtex "${tex%.*}"
pdflatex $tex
pdflatex $tex
rm "${tex%.*}.aux" "${tex%.*}.blg" "${tex%.*}.bbl" "${tex%.*}.log"


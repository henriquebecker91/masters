#!/bin/bash

name='slides_csbc'
pdflatex "$name.tex"
bibtex 	 "$name"
pdflatex "$name.tex"
# Execute a second time to create the outline
pdflatex "$name.tex"
rm *.aux *.log *.toc *.nav *.out *.snm *.bbl *.blg


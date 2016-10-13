#!/bin/sh

pdflatex ppgc-diss.tex
bibtex ppgc-diss.aux
pdflatex ppgc-diss.tex
pdflatex ppgc-diss.tex
rm *.log *.aux *.toc *.lof 


#!/bin/sh

pdflatex sbc-template.tex
bibtex sbc-template.aux
pdflatex sbc-template.tex
pdflatex sbc-template.tex
rm sbc-template.aux


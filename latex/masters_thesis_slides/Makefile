DATADIR = data
_DATA = hi_n5000-0-s731232778c5052835.csv nsds2_n50000wmin50000-0-s155213243c1070927.csv saw_n10000wmin10000-0-s985850175c608451.csv sc_a5n10000wmin10000-0-c6020768.csv ss2_wmin1000wmax1000000n10000-0-s827183242c7523725.csv
DATA = $(patsubst %,$(DATADIR)/%,$(_DATA))

all: slides_sea2016.pdf

TEXREQ = inf.eps splncs03.bst

slides_sea2016.pdf: slides_sea2016.tex instances_pya.tex $(TEXREQ) Makefile biblio.bib 
	pdflatex slides_sea2016.tex /dev/null
	bibtex slides_sea2016.aux > /dev/null
	pdflatex slides_sea2016.tex > /dev/null
	pdflatex slides_sea2016.tex > /dev/null

instances_pya.tex: instances_pya.Rnw $(DATA) Makefile
	R CMD Sweave `basename $<`

cleanaux: Makefile
	rm *.log *.aux *.toc *.lof *.out *.bbl *.lot *.blg

clean: cleanaux
	rm *.pdf


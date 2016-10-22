CHAPDIR = chapters
_CHAPTERS = instance_classes.tex
CHAPTERS = $(patsubst %,$(CHAPDIR)/%,$(_CHAPTERS))

DATADIR = data
_DATA = cutstock_knap_solvers.csv
DATA = $(patsubst %,$(DATADIR)/%,$(_DATA))

all: ppgc-diss.pdf

TEXREQ = brazilian.babel diss.def english.babel iidefs.sty iiufrgs.cls nominata.def ppgc.def

ppgc-diss.pdf: ppgc-diss.tex $(TEXREQ) $(CHAPTERS) 
	pdflatex ppgc-diss.tex > /dev/null
	bibtex ppgc-diss.aux > /dev/null
	pdflatex ppgc-diss.tex > /dev/null
	pdflatex ppgc-diss.tex > /dev/null

$(CHAPDIR)/%.tex: $(CHAPDIR)/%.Rnw $(DATA)
	cd $(CHAPDIR) && R CMD Sweave `basename $<`

cleanaux:
	rm *.log *.aux *.toc *.lof *.out *.bbl *.lot *.blg

clean: cleanaux
	rm chapters/*.tex chapters/*.aux chapters/*.pdf
	rm *.pdf


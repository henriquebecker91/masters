CHAPDIR = chapters
_CHAPTERS = introduction.tex prior_work.tex algorithms.tex instance_classes.tex experiments.tex conclusions.tex appendix.tex
CHAPTERS = $(patsubst %,$(CHAPDIR)/%,$(_CHAPTERS))

DATADIR = data
_DATA = 128_16_std_breqd_all.csv cutstock_knap_solvers.csv mtus.csv pya_ukp5.csv
DATA = $(patsubst %,$(DATADIR)/%,$(_DATA))

all: ppgc-diss.pdf

TEXREQ = brazilian.babel diss.def english.babel iidefs.sty iiufrgs.cls nominata.def ppgc.def abstract.tex resumo.tex

ppgc-diss.pdf: ppgc-diss.tex $(TEXREQ) $(CHAPTERS) Makefile biblio.bib
	pdflatex ppgc-diss.tex /dev/null
	bibtex ppgc-diss.aux > /dev/null
	pdflatex ppgc-diss.tex > /dev/null
	pdflatex ppgc-diss.tex > /dev/null

$(CHAPDIR)/%.tex: $(CHAPDIR)/%.Rnw $(DATA) Makefile
	cd $(CHAPDIR) && R CMD Sweave `basename $<`

cleanaux: Makefile
	rm *.log *.aux *.toc *.lof *.out *.bbl *.lot *.blg

clean: cleanaux
	rm chapters/*.aux chapters/*.pdf chapters/*.tex
	rm *.pdf


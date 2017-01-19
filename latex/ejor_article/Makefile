SECDIR = sections
_SECTIONS = introduction.tex prior_work.tex algorithms.tex instance_classes.tex experiments.tex conclusions.tex appendix.tex
SECTIONS = $(patsubst %,$(SECDIR)/%,$(_SECTIONS))

DATADIR = data
_DATA = 128_16_std_breqd_all.csv cutstock_knap_solvers.csv mtus.csv pya_ukp5.csv
DATA = $(patsubst %,$(DATADIR)/%,$(_DATA))

all: ejor_ukp_hb_lsb.pdf

ejor_ukp_hb_lsb.pdf: skeleton.tex $(SECTIONS) Makefile biblio.bib
	pdflatex skeleton.tex /dev/null
	bibtex skeleton.aux > /dev/null
	pdflatex skeleton.tex > /dev/null
	pdflatex skeleton.tex > /dev/null

$(SECDIR)/%.tex: $(SECDIR)/%.Rnw $(DATA) Makefile
	cd $(SECDIR) && R CMD Sweave --encoding=utf8 `basename $<`

cleanaux: Makefile
	rm -rf *.log *.aux *.toc *.lof *.out *.bbl *.lot *.blg

clean: cleanaux
	rm -rf chapters/*.aux chapters/*.pdf chapters/*.tex
	rm -rf *.pdf


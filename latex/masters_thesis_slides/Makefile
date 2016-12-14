DATADIR = data
_DATA = hi_n5000-0-s731232778c5052835.csv nsds2_n50000wmin50000-0-s155213243c1070927.csv saw_n10000wmin10000-0-s985850175c608451.csv sc_a5n10000wmin10000-0-c6020768.csv ss2_wmin1000wmax1000000n10000-0-s827183242c7523725.csv 128_16_std_breqd-n16384-s0.csv
DATA = $(patsubst %,$(DATADIR)/%,$(_DATA))

MOREDATADIR = ../../masters/thesis
_MOREDATA = 128_16_std_breqd_all.csv cutstock_knap_solvers.csv env_influence.csv mtu_impl_desktop_uc1.csv mtus.csv pya_ukp5.csv
MOREDATA = $(patsubst %,$(MOREDATADIR)/%,$(_MOREDATA))

IMAGES = $(DATA:csv=pdf) $(DATADIR)/ALL.pdf $(DATADIR)/breq_exp.pdf $(DATADIR)/env_influence.pdf $(DATADIR)/mtu1.pdf $(DATADIR)/mtu2.pdf $(DATADIR)/PP.pdf $(DATADIR)/pricing_percentage.pdf $(DATADIR)/pricing_time.pdf $(DATADIR)/Rplots.pdf $(DATADIR)/SAW.pdf $(DATADIR)/SC.pdf $(DATADIR)/SS.pdf $(DATADIR)/WCD.pdf mtu1.pdf

all: slides_sea2016.pdf

TEXREQ = inf.eps splncs03.bst

slides_sea2016.pdf: slides_sea2016.tex instances_pya.tex $(TEXREQ) Makefile biblio.bib $(IMAGES)
	pdflatex slides_sea2016.tex /dev/null
	bibtex slides_sea2016.aux > /dev/null
	pdflatex slides_sea2016.tex > /dev/null
	pdflatex slides_sea2016.tex > /dev/null

$(DATADIR)/%.pdf: $(DATADIR)/%.csv $(DATADIR)/create_images.R Makefile
	cd $(DATADIR) && R --vanilla < create_images.R

cleanaux: Makefile
	rm *.log *.aux *.toc *.lof *.out *.bbl *.lot *.blg

clean: cleanaux
	rm slides_sea2016.pdf data/*.pdf


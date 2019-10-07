all: reworking_enrichment_lit.pdf clean

clean:
	rm -f *.aux *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt
	rm -f *.fot *.cb *.cb2 .*.lb *.bcf *.run.xml *.bbl *.blg *.tex.bak

.PHONY = all clean

reworking_enrichment_lit.pdf: reworking_enrichment_lit.tex reworking_enrichment.bib
	pdflatex	reworking_enrichment_lit.tex
	biber		reworking_enrichment_lit.bcf
#	biber		reworking_enrichment_lit.bcf
	pdflatex	reworking_enrichment_lit.tex
	pdflatex	reworking_enrichment_lit.tex
	

all: barefoot_manuscript_reworking.pdf clean

clean:
	rm -f *.lof *.log *.lot *.fls *.out *.toc *.fmt
	rm -f *.fot *.cb *.cb2 .*.lb *.bcf *.run.xml *.bbl *.blg *.tex.bak
	# rm -f *.aux

remove_tex:
	rm -f barefoot_manuscript_reworking.tex *.aux

figures:
	cp ../../../figures/outputs/manuscript_figures/* figures/

.PHONY: all clean fix_biber figures remove_tex

barefoot_manuscript_reworking.pdf: barefoot_manuscript_reworking.rnw  reworking_enrichment.bib makefile response_to_reviewers.tex
	Rscript --vanilla knitit.r barefoot_manuscript_reworking.rnw barefoot_manuscript_reworking.tex
	make figures
	pdflatex	barefoot_manuscript_reworking.tex
	# pdflatex	-interaction batchmode	barefoot_manuscript_reworking.tex
	biber			barefoot_manuscript_reworking.bcf
	pdflatex	-interaction batchmode	barefoot_manuscript_reworking.tex
	pdflatex	-interaction batchmode	barefoot_manuscript_reworking.tex
	pdflatex	-interaction batchmode	response_to_reviewers.tex
	make barefoot_manuscript_reworking.docx
	make clean
	make remove_tex

barefoot_manuscript_reworking.docx:  barefoot_manuscript_reworking.tex makefile #reworking_enrichment.bib
	pandoc --reference-doc=pandoc_reference_doc.docx -o barefoot_manuscript_reworking.docx barefoot_manuscript_reworking.tex --bibliography=reworking_enrichment.bib --filter pandoc-crossref

fix_biber:
	rm -rf `biber --cache`

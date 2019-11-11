all: barefoot_manuscript_reworking.pdf clean

clean:
	rm -f *.aux *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt
	rm -f *.fot *.cb *.cb2 .*.lb *.bcf *.run.xml *.bbl *.blg *.tex.bak

figures: 
	cp ../../../figures/outputs/manuscript_figures/* figures/

.PHONY: all clean fix_biber figures

barefoot_manuscript_reworking.pdf: barefoot_manuscript_reworking.rnw  reworking_enrichment.bib makefile 
	Rscript --vanilla knitit.r barefoot_manuscript_reworking.rnw barefoot_manuscript_reworking.tex
	make figures
	pdflatex	barefoot_manuscript_reworking.tex
	biber 		barefoot_manuscript_reworking.bcf
	pdflatex	barefoot_manuscript_reworking.tex
	pdflatex	barefoot_manuscript_reworking.tex
	make clean

barefoot_manuscript_reworking.docx:  barefoot_manuscript_reworking.tex #reworking_enrichment.bib
	pandoc --reference-doc=pandoc_reference_doc.docx -o barefoot_manuscript_reworking.docx barefoot_manuscript_reworking.tex

fix_biber:
	rm -rf `biber --cache`

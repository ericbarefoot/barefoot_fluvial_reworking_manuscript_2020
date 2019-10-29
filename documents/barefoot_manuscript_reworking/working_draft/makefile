all: barefoot_manuscript_reworking.pdf clean

clean:
	rm -f *.aux *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt
	rm -f *.fot *.cb *.cb2 .*.lb *.bcf *.run.xml *.bbl *.blg *.tex.bak

.PHONY = all clean

barefoot_manuscript_reworking.pdf: barefoot_manuscript_reworking.tex  reworking_enrichment.bib
	pdflatex	barefoot_manuscript_reworking.tex
	biber 		barefoot_manuscript_reworking.bcf
	pdflatex	barefoot_manuscript_reworking.tex
	# pdflatex	barefoot_manuscript_reworking.tex
	make clean

barefoot_manuscript_reworking.docx:  barefoot_manuscript_reworking.tex #reworking_enrichment.bib
	pandoc --reference-doc=pandoc_reference_doc.docx -o barefoot_manuscript_reworking.docx barefoot_manuscript_reworking.tex

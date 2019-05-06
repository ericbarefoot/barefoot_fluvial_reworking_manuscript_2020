all: lit open

clean: 
	rm			*.bbl *.log *.blg *.aux *.xml *.bcf

lit:
	pdflatex	PROJECTNAME_lit.tex
	biber		PROJECTNAME_lit.bcf
#	biber		PROJECTNAME_lit.bcf
#	pdflatex	PROJECTNAME_lit.tex
	pdflatex	PROJECTNAME_lit.tex
	
open:
	open		headwater_stream_widths_lit.pdf
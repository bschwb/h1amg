PDF_FILES = smoothed_prol_h1.pdf smoothed_prol_hcurl.pdf

all: $(PDF_FILES)

%.pdf: %.tex
	texfot pdflatex $^ $@
	texfot pdflatex $^ $@

.PHONY: clean
clean:
	rm *.log *.aux *.pdf *.out

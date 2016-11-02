default: pdf
#default: pdf diff
################################################
# LaTeX stuff
SRC = PerrinetAdamsFriston14

edit_linux: linux_edit
linux_edit:
	texmaker $(SRC).tex &
	gedit Makefile &

edit_mac: mac_edit
mac_edit:
	mvim $(SRC).tex
	open -e Makefile

#################################################
pdf: $(SRC).pdf
#LATEXMK = latexmk -pdf
LATEXMK = latexmk --pdf  -pdflatex=pdflatex

$(SRC).pdf: $(SRC).tex
	$(LATEXMK) $(SRC)

$(SRC).bib: bib
bib: bib_citeulike
bib_citeulike:
	wget "http://www.citeulike.org/bibtex/user/LaurentPerrinet/tag/freemove" -O ../biblio/${SRC}.bib
	#wget "http://www.citeulike.org/bibtex/user/LaurentPerrinet/tag/freemove?fieldmap=posted-at:date-added&do_username_prefix=0&key_type=4" -O ../biblio/${SRC}.bib

################################################
# post-production
diff: $(SRC).tex
	#latexdiff ../13-01-02_BICY_rev0/tex/$(SRC).tex  $(SRC).tex > diff.tex
	#latexdiff ../13-11-15_BICY_rev1/$(SRC).tex  $(SRC).tex > diff.tex
	latexdiff ../14-02-06_BICY_rev2/$(SRC).tex  $(SRC).tex > diff.tex
	$(LATEXMK) diff.tex
	open diff.tex

extract_bib:
	biber --output_format=bibtex $(SRC).bcf
################################################
# cleaning macros
touch:
	touch *.tex

clean_luatex_cache:
	rm -fr ~/Library/texlive
	luaotfload-tool --update --force -vv
clean:
	rm -f diff.* *.bcf *.run.xml *.dvi *.ps *.out *.log *.aux *.bbl *.blg *.snm *.fls *.nav *.toc *.fff *.synctex.gz* *.fdb_latexmk

.PHONY:  all clean

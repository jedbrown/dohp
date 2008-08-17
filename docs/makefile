#-------------------------------------------------------
# $Id: Makefile,v 1.6 2003/09/03 19:38:20 norman Exp $
#
# This Makefile has special support for the PDF generation
# requirements for NSF proposals
# 
# This is a basic Makefile for compiling LaTeX documents
# Included is basic support for BibTeX bibliographies
# as well as support for multiple output formats.
#
# Current supported output types:
# dvi postscript html pdf
#
# Invoke the makefile in the following manners:
# make             -- Run LaTeX and BibTeX
# make dvi         -- Run LaTeX and BibTeX and invoke the DVI viewer
# make postscript  -- Make simple postscript output
# make pdf         -- Make PDF output using pdfTeX
# make finalps     -- Make better postscript output with embedded fonts
# make finalpdf    -- Make PDF output with embedded fonts, etc... 
# make index       -- Runs indexing (if you do indexing that is)
# make html        -- Makes HTML output
# make clean       -- Removes temp files and cleans up output
# make proper      -- More exhaustive clean 
#-------------------------------------------------------
#
# Name of the Main LaTeX document
# (without the .tex extenstion)
MAINDOC	= design
#
# Program Definitions
#
TEX	 = latex
BIB	 = bibtex
HTML	 = latex2html
DVIPS	 = dvips
PDFLATEX = pdflatex
PS2PDF   = ps2pdf
XDVI	 = xdvi
MAKEINDEX = makeindex
#
# Program Flags
#
TEX_FLAGS	=
BIB_FLAGS	=
HTML_FLAGS	=
PDF_FLAGS	=
PS2PDF_FLAGS    = -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true
MAKEINDEX_FLAGS = 
#DVIPS_FLAGS	= -Pcmz
DVIPS_FLAGS	= -Ppdf -G0
#
# LaTeX Source Paths
#
SRCDIR	= .
FIGDIR	= ./figures
#
# LaTeX Source Files
#
TEXSRC	= $(wildcard $(SRCDIR)/*.tex) \
	  $(wildcard $(SRCDIR)/Introduction/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_1/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_2/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_3/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_4/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_5/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_6/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_7/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_8/*.tex) \
	  $(wildcard $(SRCDIR)/Chapter_9/*.tex) 
FIGSRC	= $(wildcard $(FIGDIR)/*.ps ) \
	  $(wildcard $(FIGDIR)/*.eps) \
	  $(wildcard $(FIGDIR)/*.fig)
BIBSRC	= $(wildcard $(SRCDIR)/*.bib)
STYSRC	= $(wildcard $(SRCDIR)/*.sty) $(wildcard $(SRCDIR)/*.cls)
#
# If the Bibsrc is not empty
# then we need to generate the bbl database
# So we define it here
#
ifneq ($(strip $(BIBSRC)),)
BBLSRC	= $(MAINDOC).bbl
endif
#=====================================================#
# Display Codes (this is so we can track passes)
#=====================================================#
MOVE_TO_COL    = @echo -en "\\033[60G"
SETCOLOR_BLACK = @echo -en "\017"
SETCOLOR_RED   = @echo -en "\033[31m"
SETCOLOR_GREEN = @echo -en "\033[32m"
SETCOLOR_BLUE  = @echo -en "\033[34m"
#SETCOLOR_BLACK	= @echo -en "\\033[0;30m"
#SETCOLOR_RED 	= @echo -en "\\033[0;31m"
#SETCOLOR_GREEN 	= @echo -en "\\033[0;32m"
#SETCOLOR_BLUE	= @echo -en "\\033[0;34m"
#=====================================================#
# Standard Targets                                    #
#=====================================================#
all :	$(MAINDOC).dvi $(TEXSRC) $(FIGSRC) $(BIBSRC)

clean :	
	@rm -f *.aux *.bbl *.blg *.log

proper :
	@rm -f *.aux *.bbl *.blg *.log *.dvi \
	*.idx *.ilg *.ind *.toc *.lot *.lof 

#-----------------------------------------------------
# Portable Document Format (PDF) Output
#-----------------------------------------------------
pdf :        $(MAINDOC).pdf $(FIGSRC)
#-----------------------------------------------------
# Postscript Output
#-----------------------------------------------------
ps :         $(MAINDOC).ps $(FIGSRC)
postscript : $(MAINDOC).ps $(FIGSRC)
finalps:
	${DVIPS} ${DVIPS_FLAGS} -o ${MAINDOC}.ps ${MAINDOC}.dvi
finalpdf:
	${PS2PDF} ${PS2PDF_FLAGS} ${MAINDOC}.ps ${MAINDOC}.pdf
#-----------------------------------------------------
# Index Generation
#-----------------------------------------------------
index:	$(MAINDOC).aux $(MAINDOC).idx $(MAINDOC).ilg $(MAINDOC).ind
#-----------------------------------------------------
# HTML Output
#-----------------------------------------------------
html :       $(MAINDOC).html $(FIGSRC)
#-----------------------------------------------------
# DVI Output (with X display)
#-----------------------------------------------------
dvi :	$(MAINDOC).dvi
#	@clear
	@$(XDVI) $(MAINDOC) 1> /dev/null &
#	@clear

#=====================================================
# Compilation Rules for LaTeX and BibTeX
#=====================================================

# Note: LaTeX must be run multiple times 
# to get the proper cross referencing from
# the \ref and \cite commands
# To accomplish this we chain from tex->aux->dvi

# To generate a .aux file from a .tex file
%.aux :	%.tex
#	$(MOVE_TO_COL)
	$(SETCOLOR_GREEN)
	@echo "=========================TEX -> AUX PASS=============================="
	$(SETCOLOR_BLACK)
	@$(TEX) $(TEX_FLAGS) $(*F)
#	@rm -f $(MAINDOC).dvi

# To generate a .dvi file from a .tex file
%.dvi :	%.aux
#	$(MOVE_TO_COL)
	$(SETCOLOR_BLUE)
	@echo "=========================AUX -> DVI PASS=============================="
	$(SETCOLOR_BLACK)
	@$(TEX) $(TEX_FLAGS) $(*F)


# To generate a .idx file from a .tex file
%.ilg :	%.idx
#	$(MOVE_TO_COL)
	$(SETCOLOR_BLUE)
	@echo "=========================Indexing Pass================================"
	$(SETCOLOR_BLACK)
	@$(MAKEINDEX) $(MAKEINDEX_FLAGS) $(*F)
	@$(TEX) $(TEX_FLAGS) $(*F)

# To generate a .bbl file from a .tex file
# %.bbl : %.bib
#	@$(TEX) $(TEX_FLAGS) $(*F)
#	@$(BIB) $(BIB_FLAGS) $(*F)
#	@$(TEX) $(TEX_FLAGS) $(*F)

# To generate a .bbl file from a .tex file
%.bbl : %.tex
ifneq ($(strip $(BIBSRC)),)
#	$(MOVE_TO_COL)
	$(SETCOLOR_RED)
	@echo "===========================BIBTEX PASS================================"
	$(SETCOLOR_BLACK)
	@$(BIB) $(BIB_FLAGS) $(*F)
	$(SETCOLOR_RED)
	@echo "===========================BIBTEX/LaTeX PASS================================"
	$(SETCOLOR_BLACK)
	@$(TEX) $(TEX_FLAGS) $(*F)
endif

# To generate a .ps file from a .dvi file
%.ps : %.dvi
	@$(DVIPS) $(DVIPS_FLAGS) -o $(*F).ps $(*F).dvi

# To generate .html files from a .dvi file
%.html : %.dvi
	@$(HTML) $(HTML_FLAGS) $(*F)	

# To generate a .pdf file from a .tex file
%.pdf : %.tex
	@$(PDFLATEX) $(PDF_FLAGS) $(*F)
ifneq ($(strip $(BIBSRC)),)
	@$(BIB)      $(BIB_FLAGS) $(*F)
endif
	@$(PDFLATEX) $(PDF_FLAGS) $(*F) 
	@$(PDFLATEX) $(PDF_FLAGS) $(*F) 
#
# Dependencies
# DO NOT DELETE
$(MAINDOC).tex : $(TEXSRC) $(FIGSRC) $(STYSRC)
$(MAINDOC).aux : $(TEXSRC) $(FIGSRC) $(STYSRC)
$(MAINDOC).bbl : $(TEXSRC) $(STYSRC) $(BIBSRC)
$(MAINDOC).pdf : $(TEXSRC) $(FIGSRC) $(STYSRC) $(BBLSRC)
$(MAINDOC).dvi : $(TEXSRC) $(FIGSRC) $(STYSRC) $(BBLSRC)

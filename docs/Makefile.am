###########################################################
# Autoconf substitutions
###########################################################
LATEX=@LATEX@
DVIPS=@DVIPS@
DVIPDF=@DVIPDF@
DOXYGEN=@DOXYGEN@

#
# Rules for building dvi ps and pdf out of tex
#
SUFFIXES =.tex .dvi .ps .pdf

.tex.dvi:
	$(LATEX) $< 
	$(LATEX) $< 

.dvi.ps:
	$(DVIPS) $< -o $@

.dvi.pdf:
	$(DVIPDF) $<

#
# Install time docdir
#
docdir = @prefix@/docs

#
# Documentation targets -- apart from doxygen, the rules above should be able
# to make everything
#
doc_DATA = manual.ps manual.pdf manual.dvi doxygen

#
# This will run doxygen
#
doxygen:
	$(DOXYGEN) chromadoc

#
# Modify the install hook, so we can install directories
#
install-docDATA:
	mkdir -p $(docdir)
	cp *.ps *.dvi *.pdf $(docdir)
	cp -r ./doxygen $(docdir)/doxygen

#
# Likewise uninstall
#
uninstall-docDATA:
	cd $(docdir) ; rm -rf *.ps *.pdf *.dvi
	cd $(docdir) ; rm -rf doxygen

#
# Local clean
#
clean: 
	rm -rf doxygen
	rm -rf *.toc *.lof *.log *.aux *.ps *.dvi *.pdf

#
# Clean everything (uninstall and then local clean)
#
distclean: clean uninstall-docDATA 

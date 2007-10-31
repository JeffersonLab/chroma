#
# Recurse down into lib (to build the library)
#                   examples (build examples)
#                   includes (to install includes)
#
SUBDIRS = other_libs lib mainprogs

#
# The end user's configuration script
# 
bin_SCRIPTS = chroma-config

#
# All the other goodies
#
EXTRA_DIST = install_scripts tests docs


# xmldiff targets
.PHONY: xcheck report

xcheck:
	perl $(top_srcdir)/scripts/run_chroma_xmldiff.pl $(top_srcdir) $(top_builddir) @XMLDIFF@ $(top_builddir)/RUN

report:
	perl $(top_srcdir)/scripts/check_regres_runner.pl $(top_srcdir) $(top_builddir) @XMLDIFF@
	xsltproc $(top_srcdir)/scripts/regres_report_to_html.xsl $(top_builddir)/regres_report.xml > regres_report.html



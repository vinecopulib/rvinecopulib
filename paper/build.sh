#!/bin/sh
# Build the manuscript. R ships jss.cls/jss.bst under its own texmf tree.
# The real jsslogo ships with the JSS author kit.  For local drafts we make a
# 1pt placeholder so the document compiles; replace it before submitting.
if [ ! -f jsslogo.pdf ]; then
  echo "note: jsslogo.pdf missing, generating a placeholder (replace before submission)"
  printf '%s\n' '\documentclass{article}' \
    '\usepackage[paperwidth=1pt,paperheight=1pt,margin=0pt]{geometry}' \
    '\pagestyle{empty}\begin{document}\hbox{}\end{document}' > _logo.tex
  pdflatex -interaction=batchmode _logo.tex >/dev/null 2>&1 && mv _logo.pdf jsslogo.pdf
  rm -f _logo.tex _logo.aux _logo.log
fi

R_TEXMF=$(Rscript -e 'cat(file.path(R.home("share"), "texmf"))' 2>/dev/null)
export TEXINPUTS="${R_TEXMF}/tex/latex:${TEXINPUTS}"
export BSTINPUTS="${R_TEXMF}/bibtex/bst:${BSTINPUTS}"
pdflatex -interaction=nonstopmode rvinecopulib.tex
bibtex rvinecopulib
pdflatex -interaction=nonstopmode rvinecopulib.tex
pdflatex -interaction=nonstopmode rvinecopulib.tex

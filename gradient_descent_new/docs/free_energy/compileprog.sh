#! /bin/bash

# compiles tex file with bibiography, suppresses warnings but exits
# if an error occurs

mkdir -p intfiles # if no subdirectory to store intermediate files, then
#                   make the subdirectory

pdflatex -halt-on-error --output-directory=intfiles $1.tex\
| grep -a3 ^! # grep -a3 ^! checks if there is a ! three lines up from error
greprc=$? # if grep finds something (i.e. error occurs, then return value
#           is currently at 0, if it doesn't find anything then return
#           value is 1 (unless something weird happens and return is 2)

if [[ $greprc -eq 0 ]] ; then
    exit 1
fi
echo First compilation of pdflatex successful
bibtex intfiles/$1.aux
pdflatex -halt-on-error --output-directory=intfiles $1.tex\
| grep -a3 ^!
greprc=$?
if [[ $greprc -eq 0 ]] ; then
    exit 1
fi
echo Second compilation of pdflatex successful
pdflatex -halt-on-error --output-directory=intfiles $1.tex\
| grep -a3 ^!
greprc=$?
if [[ $greprc -eq 0 ]] ; then
    exit 1
fi
echo Final compilation of pdflatex successful
# move the pdf to the cwd

mv intfiles/$1.pdf $1.pdf
#!/bin/bash
echo 'Calling latex'
latex pypedal > /dev/null
echo 'Calling bibtex'
bibtex pypedal > /dev/null
echo 'Calling latex'
latex pypedal > /dev/null
echo 'Calling latex'
latex pypedal > /dev/null
echo 'Calling makeindex'
makeindex pypedal.idx > /dev/null
echo 'Calling makeindex'
makeindex -o pypedal.fnd pypedal.fdx > /dev/null
echo 'Calling latex'
latex pypedal > /dev/null
echo 'Calling dvips'
dvips -t letter -o pypedal.ps pypedal > /dev/null
echo 'Calling ps2pdf'
ps2pdf pypedal.ps pypedal.pdf > /dev/null
echo 'Calling latex2html'
tools/mkhowto --html --about html/stdabout.dat --iconserver ../icons --favicon ../icons/pyfav.png --address "See <i><a href=\"about.html\">About this document...</a></i> for information on suggesting changes." --up-link ../index.html --up-title "Python Documentation Index" --global-module-index "../modindex.html" --dvips-safe pypedal.tex

# Now we have to use perl to change the folder in which icons are located
# because the Apache installation at Sourceforge redorects /icons to a
# system folder. If we don;t do this, the navication buttons in the HTML
# version are broken.
echo 'Calling perl to fix icon paths'
perl -pi -e 's/..\/icons\//manicons\//g' /home/jcole/pypedal/PyPedal/doc/pypedal/*.html

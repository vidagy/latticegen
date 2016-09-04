#!/bin/bash
filename=$1
echo $filename

plotFN="plot_${filename}"

gnuplot <<-EOP2
set terminal epslatex standalone color solid size 5,4
set size ratio 1.0
set angles degrees
set key t l  Left

set xlabel "\$ x \$"
set ylabel "\$ y \$"

set xrange [0:1.0]
set yrange [0:1.0]

#set format x "% .1f"
#set format y "% .2f"

set grid

set output "${plotFN}.tex"

plot "${filename}" u 1:2 w p pt 7 lc rgb "red" notitle

EOP2

latex ${plotFN}.tex
dvipdf ${plotFN}.dvi
rm ${plotFN}.log ${plotFN}.aux ${plotFN}.tex ${plotFN}-inc.eps ${plotFN}.dvi fit.log





 

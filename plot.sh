#!/bin/bash

N="data/tmp"
PLOTNM="$N.eps"
DATA="$N.txt"
SIZE=(`stty size`)
gnuplot -e "nm='$PLOTNM'; data='$DATA'; size1='${SIZE[1]}'; size2='${SIZE[0]}'" gnuplot_params.p
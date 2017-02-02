# output ASCII plot directly to terminal
set term dumb size1,size2 
#if (!exists("nm")) nm='plot.eps'
#set output nm
#set xrange [0.95:1]
plot data using 1:4
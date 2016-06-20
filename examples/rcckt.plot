#set xrange [0:1e-6]
plot "rcckt.tran" using 1:2 with lines, \
 "rcckt.tran" using 1:3 with lines, \
 "spice.dat" using 1:2 with lines, \
 "spice.dat" using 1:3 

pause -1 "Hit return to exit gnuplot"

#set xrange [0:1e-6]
plot "example4.tran" using 1:2 with lines, \
 "example4.tran" using 1:3 with lines, \
 "example4_spice.dat" using 1:3 with lines

pause -1 "Hit return to exit gnuplot"

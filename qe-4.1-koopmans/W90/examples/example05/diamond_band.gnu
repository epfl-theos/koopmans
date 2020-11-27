set style data dots
set nokey
set xrange [0: 6.38491]
set yrange [ -7.43859 : 20.39778]
set arrow from  1.68570,  -7.43859 to  1.68570,  20.39778 nohead
set arrow from  3.63217,  -7.43859 to  3.63217,  20.39778 nohead
set arrow from  4.32036,  -7.43859 to  4.32036,  20.39778 nohead
set xtics (" L "  0.00000," G "  1.68570," X "  3.63217," K "  4.32036," G "  6.38491)
 plot "diamond_band.dat"

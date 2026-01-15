set term pdfcairo dashed enhanced size 4, 2
set datafile separator " "

set grid
set xrange [1:8]
#set yrange [1e-9:0.2]
set xtics
#set xtics add (5, 25, 50)
set ytics
#set logscale x
#set logscale y
#set format y "10^{%L}"
#set ytics 0.002
set xlabel "Number of cores"
set ylabel "Normalised execution time (speed-up)"
set key right bottom outside;
set rmargin 21
set key spacing 1.2

set label "Dashed lines are speed-up \n for looping only" at graph 1.05,0.9 font ",10"


# Data lines for each p-order (with titles)
set style line 1 lc rgb "black"        pt 7 ps 0.5 lw 1
set style line 2 lc rgb "slategrey"    pt 7 ps 0.5 lw 1
set style line 3 lc rgb "pink"         pt 7 ps 0.5 lw 1
set style line 4 lc rgb "red"          pt 7 ps 0.5 lw 1
set style line 5 lc rgb "blue"         pt 7 ps 0.5 lw 1
set style line 6 lc rgb "violet"       pt 7 ps 0.5 lw 1

set style line 11 lc rgb "black"        pt 7 ps 0.5 lw 1 dt 2
set style line 22 lc rgb "slategrey"    pt 7 ps 0.5 lw 1 dt 2
set style line 33 lc rgb "pink"         pt 7 ps 0.5 lw 1 dt 2
set style line 44 lc rgb "red"          pt 7 ps 0.5 lw 1 dt 2
set style line 55 lc rgb "blue"         pt 7 ps 0.5 lw 1 dt 2
set style line 66 lc rgb "violet"       pt 7 ps 0.5 lw 1 dt 2

set style line 10 lc rgb "grey" dt 3 lw 2

set output "strongScaling_hex_2D.pdf
set title "2D hexahedral mesh"

plot \
     x w l ls 10 title "Ideal Speed-up", \
    "hex-2D.mesh-1.summary.txt" u 1:4 w lp ls 1  title "400 cells", \
    "hex-2D.mesh-2.summary.txt" u 1:4 w lp ls 2  title "1 600 cells", \
    "hex-2D.mesh-3.summary.txt" u 1:4 w lp ls 3  title "6 400 cells", \
    "hex-2D.mesh-4.summary.txt" u 1:4 w lp ls 4  title "25 600 cells", \
    "hex-2D.mesh-5.summary.txt" u 1:4 w lp ls 5  title "102 400 cells", \
    "hex-2D.mesh-6.summary.txt" u 1:4 w lp ls 6  title "409 600 cells",\
    "hex-2D.mesh-1.summary.txt" u 1:6 w lp ls 11  notitle,\
    "hex-2D.mesh-2.summary.txt" u 1:6 w lp ls 22  notitle,\
    "hex-2D.mesh-3.summary.txt" u 1:6 w lp ls 33  notitle,\
    "hex-2D.mesh-4.summary.txt" u 1:6 w lp ls 44  notitle,\
    "hex-2D.mesh-5.summary.txt" u 1:6 w lp ls 55  notitle,\
    "hex-2D.mesh-6.summary.txt" u 1:6 w lp ls 66  notitle


set output "strongScaling_hex_3D.pdf
set title "3D hexahedral mesh"

plot \
     x w l ls 10 title "Ideal Speed-up", \
    "hex-3D.mesh-1.summary.txt" u 1:4 w lp ls 1  title "1 000 cells", \
    "hex-3D.mesh-2.summary.txt" u 1:4 w lp ls 2  title "4 096 cells", \
    "hex-3D.mesh-3.summary.txt" u 1:4 w lp ls 3  title "15 625 cells", \
    "hex-3D.mesh-4.summary.txt" u 1:4 w lp ls 4  title "64 000 cells", \
    "hex-3D.mesh-5.summary.txt" u 1:4 w lp ls 5  title "262 144 cells", \
    "hex-3D.mesh-6.summary.txt" u 1:4 w lp ls 6  title "1 000 000 cells",\
    "hex-3D.mesh-1.summary.txt" u 1:6 w lp ls 11  notitle,\
    "hex-3D.mesh-2.summary.txt" u 1:6 w lp ls 22  notitle,\
    "hex-3D.mesh-3.summary.txt" u 1:6 w lp ls 33  notitle,\
    "hex-3D.mesh-4.summary.txt" u 1:6 w lp ls 44  notitle,\
    "hex-3D.mesh-5.summary.txt" u 1:6 w lp ls 55  notitle,\
    "hex-3D.mesh-6.summary.txt" u 1:6 w lp ls 66  notitle


set yrange [0:8]
set output "strongScaling_hex_3D_limitedY.pdf
set title "3D hexahedral mesh"

plot \
     x w l ls 10 title "Ideal Speed-up", \
    "hex-3D.mesh-1.summary.txt" u 1:4 w lp ls 1  title "1 000 cells", \
    "hex-3D.mesh-2.summary.txt" u 1:4 w lp ls 2  title "4 096 cells", \
    "hex-3D.mesh-3.summary.txt" u 1:4 w lp ls 3  title "15 625 cells", \
    "hex-3D.mesh-4.summary.txt" u 1:4 w lp ls 4  title "64 000 cells", \
    "hex-3D.mesh-5.summary.txt" u 1:4 w lp ls 5  title "262 144 cells", \
    "hex-3D.mesh-6.summary.txt" u 1:4 w lp ls 6  title "1 000 000 cells",\
    "hex-3D.mesh-1.summary.txt" u 1:6 w lp ls 11  notitle,\
    "hex-3D.mesh-2.summary.txt" u 1:6 w lp ls 22  notitle,\
    "hex-3D.mesh-3.summary.txt" u 1:6 w lp ls 33  notitle,\
    "hex-3D.mesh-4.summary.txt" u 1:6 w lp ls 44  notitle,\
    "hex-3D.mesh-5.summary.txt" u 1:6 w lp ls 55  notitle,\
    "hex-3D.mesh-6.summary.txt" u 1:6 w lp ls 66  notitle

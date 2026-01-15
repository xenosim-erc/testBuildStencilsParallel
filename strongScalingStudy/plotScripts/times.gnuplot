set term pdfcairo dashed enhanced size 4, 2
set datafile separator " "

set grid
#set xrange [1:8]
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

set output "timings_hex_2D.pdf
set title "2D hexahedral mesh"


set style data histograms
set style histogram clustered gap 1
set style fill solid 0.8 border -1
set boxwidth 0.9

set xtics rotate by -45
set grid ytics

set xlabel "Cores"
set ylabel "Time (seconds)"
set title "Time vs Looping-Time by Number of Cores"
set xtics 1

# Use Cores as x labels
set xtics nomirror

set label 1 "400 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_400.pdf

plot \
    'hex-2D.mesh-1.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-2D.mesh-1.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "1 600 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_1600.pdf
plot \
    'hex-2D.mesh-2.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-2D.mesh-2.summary.txt' using 5:xtic(1) title "Looping-Time

set label 1 "6 400 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_6400.pdf
plot \
     'hex-2D.mesh-3.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-2D.mesh-3.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "25600 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_25600.pdf
plot \
     'hex-2D.mesh-4.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-2D.mesh-4.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "102 400 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_102400.pdf
plot \
     'hex-2D.mesh-5.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-2D.mesh-5.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "409 600 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_2D_409600.pdf
plot \
     'hex-2D.mesh-6.summary.txt' using 3:xtic(1) title "Total time", \
     'hex-2D.mesh-6.summary.txt' using 5:xtic(1) title "Looping-Time"


set label 1 "1000 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_1000.pdf

plot \
    'hex-3D.mesh-1.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-3D.mesh-1.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "4 096 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_4096.pdf
plot \
    'hex-3D.mesh-2.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-3D.mesh-2.summary.txt' using 5:xtic(1) title "Looping-Time

set label 1 "15 625 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_15625.pdf
plot \
    'hex-3D.mesh-3.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-3D.mesh-3.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "64 000 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_64000.pdf
plot \
    'hex-3D.mesh-4.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-3D.mesh-4.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "262 144 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_262144.pdf
plot \
    'hex-3D.mesh-5.summary.txt' using 3:xtic(1) title "Total time", \
    'hex-3D.mesh-5.summary.txt' using 5:xtic(1) title "Looping-Time"

set label 1 "1 000 000 cells" at graph 1.05,0.9 font ",10"
set output "timings_hex_3D_1000000.pdf
plot \
     'hex-3D.mesh-6.summary.txt' using 3:xtic(1) title "Total time", \
     'hex-3D.mesh-6.summary.txt' using 5:xtic(1) title "Looping-Time"

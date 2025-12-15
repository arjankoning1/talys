set xlabel "Incident Photon Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 90}Zr({/=22{/Symbol g}},x{/=25{n}})^{/=14 }" font ",25"
#
set pointsize 1.2
set xtics font "Times-Roman, 18" 
set ytics font "Times-Roman, 18"
#
# Setting of the axis format  
#
set format x "%G"
set format y "%G"
#
#  Setting of the legend location  
#
set key font ",16"
set key top right 
#set grid x
#set grid y
# show grid
# 
set term postscript enhanced font "Times-Roman,16" color
set output "g-Zr090-xs.eps"
#
set style line 10 linetype 1 linewidth 3 linecolor rgb "black" 
#
#
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 1 linetype 4  linewidth 6 linecolor rgb "blue"
set style line 1 linetype 4  linewidth 2 linecolor rgb "blue"
# set style line 1 linetype 4  linewidth 2 linecolor rgb "blue" dashtype '-..'
#
#######################################################################################
  set style line 2 linetype 3  linewidth 2 linecolor rgb "#C71585"
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 3 linetype 8  linewidth 6 linecolor rgb "#C71585"
set style line 3 linetype 8  linewidth 2 linecolor rgb "#C71585"
# set style line 3 linetype 8  linewidth 2 linecolor rgb "#C71585" dashtype '-..'
#
#######################################################################################
set style line 11 pointtype 1 linecolor rgb "blue "
set style line 12 pointtype 2 linecolor rgb "#00008B "
plot [0.:30][2.84992013E-04:301.760986] \
"/Users/koning/talys/samples/g-Zr090-xs/org/nprod.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS" w lines linestyle 10, \
"/Users/koning/libraries/g/Zr090/endfb8.1/tables/xs/g-Zr090-MT201.endfb8.1" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "ENDF/B-VIII.1" w lines linestyle 1, \
"/Users/koning/libraries/g/Zr090/iaea.pd/tables/xs/g-Zr090-MT201.iaea.pd" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "IAEA" w lines linestyle 2, \
"/Users/koning/libraries/g/Zr090/tendl.2023/tables/xs/g-Zr090-MT201.tendl.2023" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TENDL-2023" w lines linestyle 3, \
"/Users/koning/libraries/g/Zr090/exfor/xs/201/g-Zr090-MT201-Lepretre-L0027011.1971" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/g/Zr090/exfor/xs/201/g-Zr090-MT201-Lepretre-L0027011.1971" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Lepretre(1971) " linestyle 11, \
"/Users/koning/libraries/g/Zr090/exfor/xs/201/g-Zr090-MT201-Berman-L0011005.1967" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 12, \
"/Users/koning/libraries/g/Zr090/exfor/xs/201/g-Zr090-MT201-Berman-L0011005.1967" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Berman(1967) " linestyle 12, \

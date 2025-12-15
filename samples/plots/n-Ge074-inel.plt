set xlabel "Incident Neutron Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 74}Ge(n,{/=25{n'}})^{/=14 }" font ",25"
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
set output "n-Ge074-inel.eps"
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
# set style line 1 linetype 2  linewidth 6 linecolor rgb "red"
set style line 1 linetype 2  linewidth 2 linecolor rgb "red"
# set style line 1 linetype 2  linewidth 2 linecolor rgb "red" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 2 linetype 3  linewidth 6 linecolor rgb "#32CD32"
set style line 2 linetype 3  linewidth 2 linecolor rgb "#32CD32"
# set style line 2 linetype 3  linewidth 2 linecolor rgb "#32CD32" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 3 linetype 4  linewidth 6 linecolor rgb "blue"
set style line 3 linetype 4  linewidth 2 linecolor rgb "blue"
# set style line 3 linetype 4  linewidth 2 linecolor rgb "blue" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 4 linetype 3  linewidth 6 linecolor rgb "#2E8B57"
set style line 4 linetype 3  linewidth 2 linecolor rgb "#2E8B57"
# set style line 4 linetype 3  linewidth 2 linecolor rgb "#2E8B57" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 5 linetype 8  linewidth 6 linecolor rgb "#C71585"
set style line 5 linetype 8  linewidth 2 linecolor rgb "#C71585"
# set style line 5 linetype 8  linewidth 2 linecolor rgb "#C71585" dashtype '-..'
#
#######################################################################################
set style line 11 pointtype 1 linecolor rgb "red "
set style line 12 pointtype 2 linecolor rgb "red "
set style line 13 pointtype 3 linecolor rgb "red "
set style line 14 pointtype 4 linecolor rgb "red "
set style line 15 pointtype 5 linecolor rgb "red "
set style line 16 pointtype 6 linecolor rgb "orange "
set style line 17 pointtype 7 linecolor rgb "magenta "
set style line 18 pointtype 8 linecolor rgb "green "
set style line 19 pointtype 9 linecolor rgb "blue "
set style line 20 pointtype 10 linecolor rgb "blue "
plot [0:20][0.:1500.0] \
"/Users/koning/talys/samples/n-Ge074-vib/org/nn.L01" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "Spherical OMP" w lines linestyle 10, \
"/Users/koning/libraries/n/Ge074/exfor/xs/051/n-Ge074-MT051-Konobeevskii-40196006.1972" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/n/Ge074/exfor/xs/051/n-Ge074-MT051-Konobeevskii-40196006.1972" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) w points t "Konobeevskii(1972) " linestyle 11, \
"/Users/koning/libraries/n/Ge074/exfor/xs/051/n-Ge074-MT051-Chung-10045069.1970" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 12, \
"/Users/koning/libraries/n/Ge074/exfor/xs/051/n-Ge074-MT051-Chung-10045069.1970" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Chung(1970) " linestyle 12, \

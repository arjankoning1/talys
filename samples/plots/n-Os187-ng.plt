set xlabel "Incident Neutron Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 187}Os(n,{/=22{/Symbol g}})^{/=14 188}Os" font ",25"
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
set output "n-Os187-ng.eps"
#
set logscale xy
set format x "10^{%L}"
set format y "10^{%L"
#
set logscale
set format x "%G"
set format y "10^{%L}"
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
plot [1.e-3:1][100.:10000.0] \
"/Users/koning/talys/samples/n-Os187-astro-ng/org/xs000000.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS" w lines linestyle 10, \
"/Users/koning/libraries/n/Os187/jeff4.0/tables/xs/n-Os187-MT102.jeff4.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JEFF-3.3" w lines linestyle 1, \
"/Users/koning/libraries/n/Os187/jendl5.0/tables/xs/n-Os187-MT102.jendl5.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JENDL-5.0" w lines linestyle 2, \
"/Users/koning/libraries/n/Os187/endfb8.1/tables/xs/n-Os187-MT102.endfb8.1" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "ENDF/B-VIII.1" w lines linestyle 3, \
"/Users/koning/libraries/n/Os187/tendl.2023/tables/xs/n-Os187-MT102.tendl.2023" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TENDL-2023" w lines linestyle 5, \
"/Users/koning/libraries/n/Os187/exfor/xs/102/n-Os187-MT102-Winters-10882003.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/n/Os187/exfor/xs/102/n-Os187-MT102-Winters-10882003.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Winters(1980) " linestyle 11, \

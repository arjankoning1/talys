set xlabel "Neutron number" font ",18"
set ylabel "P" font ",18" offset 0
set bmargin 5
set title "^{/=14 239}Pu(n_{th},{/=25{f}})^{/=14 } P(nu)" font ",25"
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
set output "n-Pu239-FY-Pnu.eps"
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
# set style line 4 linetype 7  linewidth 6 linecolor rgb "#B0C4DE"
set style line 4 linetype 7  linewidth 2 linecolor rgb "#B0C4DE"
# set style line 4 linetype 7  linewidth 2 linecolor rgb "#B0C4DE" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 5 linetype 3  linewidth 6 linecolor rgb "#2E8B57"
set style line 5 linetype 3  linewidth 2 linecolor rgb "#2E8B57"
# set style line 5 linetype 3  linewidth 2 linecolor rgb "#2E8B57" dashtype '-..'
#
#######################################################################################
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 6 linetype 8  linewidth 6 linecolor rgb "#C71585"
set style line 6 linetype 8  linewidth 2 linecolor rgb "#C71585"
# set style line 6 linetype 8  linewidth 2 linecolor rgb "#C71585" dashtype '-..'
#
#######################################################################################
set style line 11 pointtype 1 linecolor rgb "#DC143C "
set style line 12 pointtype 2 linecolor rgb "#DC143C "
set style line 13 pointtype 3 linecolor rgb "#DC143C "
set style line 14 pointtype 4 linecolor rgb "#DC143C "
set style line 15 pointtype 5 linecolor rgb "#DC143C "
set style line 16 pointtype 6 linecolor rgb "red "
set style line 17 pointtype 7 linecolor rgb "#FF8C00 "
set style line 18 pointtype 8 linecolor rgb "orange "
set style line 19 pointtype 9 linecolor rgb "magenta "
set style line 20 pointtype 10 linecolor rgb "magenta "
set style line 21 pointtype 11 linecolor rgb "magenta "
set style line 22 pointtype 12 linecolor rgb "cyan "
set style line 23 pointtype 13 linecolor rgb "cyan "
set style line 24 pointtype 14 linecolor rgb "cyan "
set style line 25 pointtype 15 linecolor rgb "cyan "
set style line 26 pointtype 16 linecolor rgb "green "
set style line 27 pointtype 17 linecolor rgb "green "
set style line 28 pointtype 18 linecolor rgb "blue "
set style line 29 pointtype 19 linecolor rgb "blue "
set style line 30 pointtype 20 linecolor rgb "blue "
plot [0.:8.0][0.:0.4] \
"/Users/koning/talys/samples/n-Pu239-fy-gef/org/Pnun1.00E-06.fis" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "GEF" w boxes linestyle 10, \
"/Users/koning/talys/samples/n-Pu239-fy-hf3d/org/Pnun1.00E-06.fis" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "HF3D" w boxes linestyle 1, \

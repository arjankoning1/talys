set xlabel "Irradiation + cooling time [h]" font ",18"
set ylabel "Activity [GBq]" font ",18" offset 0
set bmargin 5
set title "^{/=14 100}Mo(p,{/=25{2n}})^{/=14 99m}Tc radioactive yield     " font ",25"  
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
set term postscript enhanced font "Times-Roman,16" color
set output "p-Mo100-yield.eps"
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
# set style line 1 linetype 3  linewidth 4 linecolor rgb "green"
set style line 1 linetype 3  linewidth 2 linecolor rgb "green"
# set style line 1 linetype 3  linewidth 2 linecolor rgb "green" dashtype '-..'
#
#######################################################################################
  set style line 2 linetype 3  linewidth 2 linecolor rgb "cyan"
#
######################################################################################
#
# If you want to have only solid lines then uncomment line below !!!
#  
#  a/ tick solid lines, b/ thinner solid lines, c/ dashed lines (acceptable for papers!) 
#
# set style line 3 linetype 8  linewidth 4 linecolor rgb "orange"
set style line 3 linetype 8  linewidth 2 linecolor rgb "orange"
# set style line 3 linetype 8  linewidth 2 linecolor rgb "orange" dashtype '-..'
#
#######################################################################################
  set style line 11 pointtype 1 linecolor rgb "red "
  set style line 12 pointtype 2 linecolor rgb "#FF4500 "
  set style line 13 pointtype 3 linecolor rgb "#FF4500 "
  set style line 14 pointtype 4 linecolor rgb "#FF4500 "
  set style line 15 pointtype 5 linecolor rgb "#FF4500 "
  set style line 16 pointtype 6 linecolor rgb "#FF4500 "
  set style line 17 pointtype 7 linecolor rgb "#FF4500 "
  set style line 18 pointtype 8 linecolor rgb "magenta "

plot [0:48][0:1000] \
"/Users/koning/talys/samples/p-Mo100-medical/org/Y043099.L02" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS" w lines linestyle 10, \

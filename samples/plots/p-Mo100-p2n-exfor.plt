set xlabel "Incident Proton Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 100}Mo(p,{/=25{2n}})^{/=14 99m}Tc     " font ",25"  
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
set output "p-Mo100-p2n-exfor.eps"
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

plot [5:30][0:500] \
"/Users/koning/talys/samples/p-Mo100-medical/org/rp043099.L02" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS" w lines linestyle 10, \
"/Users/koning/libraries/p/Mo100/jendl5.0/tables/residual/p-Mo100-rp043099m.jendl5.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JENDL-5.0" w lines linestyle 1, \
"/Users/koning/libraries/p/Mo100/iaea.2024/tables/residual/p-Mo100-rp043099m.iaea.2024" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "IAEA" w lines linestyle 2, \
"/Users/koning/libraries/p/Mo100/tendl.2023/tables/residual/p-Mo100-rp043099m.tendl.2023" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TENDL-2023" w lines linestyle 3, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Lamere-C2413034.2019" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Lamere-C2413034.2019" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 11, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156005.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 12, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156005.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 12, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156006.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 13, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156006.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 13, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156007.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 14, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Gagnon-C2156007.2011" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 14, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Takacs-D4322002.2015" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 15, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Takacs-D4322002.2015" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 15, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Takacs-D4322004.2015" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 16, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Takacs-D4322004.2015" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 16, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Tarkanyi-D4264010.2012" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 17, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Tarkanyi-D4264010.2012" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 17, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Lagunassolar-C0963002.1996" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 18, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Lagunassolar-C0963002.1996" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 18, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Levkovski-A0510249.1991" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 19, \
"/Users/koning/libraries/p/Mo100/exfor/residual/043099m/p-Mo100-rp043099m-Levkovski-A0510249.1991" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "EXFOR " linestyle 19, \


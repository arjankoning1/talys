set xlabel "Incident Neutron Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 232}Th(n,{/=25{f}})^{/=14 }" font ",25"
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
set output "n-Th232-fis-wkb.eps"
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
plot [0.:20.0000000][0.:600.0] \
"/Users/koning/talys/samples/n-Th232-fis-wkb/org/fission.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS" w lines linestyle 10, \
"/Users/koning/libraries/n/Th232/jeff4.0/tables/xs/n-Th232-MT018.jeff4.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JEFF-3.3" w lines linestyle 1, \
"/Users/koning/libraries/n/Th232/jendl5.0/tables/xs/n-Th232-MT018.jendl5.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JENDL-5.0" w lines linestyle 2, \
"/Users/koning/libraries/n/Th232/endfb8.1/tables/xs/n-Th232-MT018.endfb8.1" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "ENDF/B-VIII.1" w lines linestyle 3, \
"/Users/koning/libraries/n/Th232/irdff2.0/tables/xs/n-Th232-MT018.irdff2.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "IRDFF-2.0" w lines linestyle 4, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Jain-31424002.1997" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 18, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Jain-31424002.1997" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Jain(1997) " linestyle 18, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Garlea-31459025.1992" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 19, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Garlea-31459025.1992" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Garlea(1992) " linestyle 19, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Anand-30975002.1986" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 22, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Anand-30975002.1986" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Anand(1986) " linestyle 22, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Lisowski-14176006.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 23, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Lisowski-14176006.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Lisowski(1988) " linestyle 23, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Manabe-222820022.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 24, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Manabe-222820022.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Manabe(1988) " linestyle 24, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Meadows-131340032.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 25, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Meadows-131340032.1988" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Meadows(1988) " linestyle 25, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Meadows-108430032.1979" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 26, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Meadows-108430032.1979" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Meadows(1983) " linestyle 26, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Blons-21656003.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 29, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Blons-21656003.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Blons(1980) " linestyle 29, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Casanova-20953002.1973" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 30, \
"/Users/koning/libraries/n/Th232/exfor/xs/018/n-Th232-MT018-Casanova-20953002.1973" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Casanova(1973) " linestyle 30, \

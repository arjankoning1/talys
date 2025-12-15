set xlabel "Incident Alpha Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 165}Ho({/=22{/Symbol a}},{/=25{n}})^{/=14 168}Tm" font ",25"
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
set output "a-Ho165-omp.eps"
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
# set style line 1 linetype 1  linewidth 6 linecolor rgb "black"
set style line 1 linetype 1  linewidth 2 linecolor rgb "black"
# set style line 1 linetype 1  linewidth 2 linecolor rgb "black" dashtype '-..'
#
#######################################################################################
set style line 11 pointtype 1 linecolor rgb "#DC143C "
set style line 12 pointtype 2 linecolor rgb "red "
set style line 13 pointtype 3 linecolor rgb "orange "
set style line 14 pointtype 4 linecolor rgb "magenta "
set style line 15 pointtype 5 linecolor rgb "magenta "
set style line 16 pointtype 6 linecolor rgb "magenta "
set style line 17 pointtype 7 linecolor rgb "blue "
set style line 18 pointtype 8 linecolor rgb "#00008B "
set style line 19 pointtype 9 linecolor rgb "#00008B "
set style line 20 pointtype 10 linecolor rgb "black "
plot [1.:50][1.95032808E-05:100.] \
"/Users/koning/talys/samples/a-Ho165-omp1/org/rp069168.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "alphaomp 1" w lines linestyle 1, \
"/Users/koning/talys/samples/a-Ho165-omp2/org/rp069168.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "alphaomp 2" w lines linestyle 2, \
"/Users/koning/talys/samples/a-Ho165-omp5/org/rp069168.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "alphaomp 5" w lines linestyle 3, \
"/Users/koning/talys/samples/a-Ho165-omp6/org/rp069168.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "alphaomp 6" w lines linestyle 4, \
"/Users/koning/libraries/a/Ho165/tendl.2023/tables/xs/a-Ho165-MT004.tendl.2023" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) t "TENDL-2023" w lines linestyle 5, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Glorius-C2119003.2014" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Glorius-C2119003.2014" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Glorius(2014) " linestyle 11, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Tarkanyi-D4275002.2010" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 12, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Tarkanyi-D4275002.2010" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Tarkanyi(2010) " linestyle 12, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Gadkari-D0068002.1997" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 13, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Gadkari-D0068002.1997" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Gadkari(1997) " linestyle 13, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Mukherjee-D4225004.1991" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 14, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Mukherjee-D4225004.1991" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Mukherjee(1991) " linestyle 14, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Singh-D4195002.1995" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 15, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Singh-D4195002.1995" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Singh(1995) " linestyle 15, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Singh-O1085002.1992" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 16, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Singh-O1085002.1992" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Singh(1992) " linestyle 16, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Homma-D4222012.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 17, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Homma-D4222012.1980" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Homma(1980) " linestyle 17, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-MartinJr-R0052004.1966" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 18, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-MartinJr-R0052004.1966" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "MartinJr(1966) " linestyle 18, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Sau-D0600002.1968" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 19, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Sau-D0600002.1968" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Sau(1968) " linestyle 19, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Wilkinson-D4215004.1949" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 20, \
"/Users/koning/libraries/a/Ho165/exfor/xs/004/a-Ho165-MT004-Wilkinson-D4215004.1949" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Wilkinson(1949) " linestyle 20, \

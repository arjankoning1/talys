set xlabel "Incident Proton Energy [MeV]" font ",18"
set ylabel "Cross section [mb]" font ",18" offset 0
set bmargin 5
set title "^{/=14 123}Te(p,{/=25{n}})^{/=14 123}I     " font ",25"  
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
set output "p-Te123-pn-fit.eps"
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
  set style line 11 pointtype 1 linecolor rgb "magenta "
  set style line 12 pointtype 2 linecolor rgb "#00FF00 "
  set style line 13 pointtype 3 linecolor rgb "#32CD32 "
plot [1.:30][1.41174402E-04:1000] \
"/Users/koning/talys/samples/p-Te123-pn-fit/org/rp053123.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS fitted" w lines linestyle 10, \
"/Users/koning/libraries/p/Te123/jendl5.0/tables/xs/p-Te123-MT004.jendl5.0" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "JENDL-5.0" w lines linestyle 1, \
"/Users/koning/libraries/p/Te123/iaea.2024/tables/residual/p-Te123-rp053123.iaea.2024" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "IAEA" w lines linestyle 2, \
"/Users/koning/talys/samples/p-Te123-pn-global/org/rp053123.tot" u ( $1 > 0. ? $1 : 1.e-10 ):( $2 > 0. ? $2 : 1.e-10 ) t "TALYS global" w lines linestyle 3, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Mahunka-D4142002.1996" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 11, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Mahunka-D4142002.1996" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Mahunka(1996) " linestyle 11, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Scholten-A04730021.1989" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 12, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Scholten-A04730021.1989" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Scholten(1989) " linestyle 12, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Barrall-D0093002.1981" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ):2:4 w xyerrorbars notitle linestyle 13, \
"/Users/koning/libraries/p/Te123/exfor/xs/004/p-Te123-MT004-Barrall-D0093002.1981" u ( $1 > 0. ? $1 : 1.e-10 ):( $3 > 0. ? $3 : 1.e-10 ) w points t "Barrall(1981) " linestyle 13, \

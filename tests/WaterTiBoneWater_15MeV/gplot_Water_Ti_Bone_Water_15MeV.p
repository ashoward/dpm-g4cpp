set term postscript eps enhanced color dl 1.1 "Times" 20
######## Output File Name
#
set output "fig_DDose_Water_Ti_Bone_Water_15MeV.eps"
######## Line types
#
set title "E_0 = 15.0 [MeV], Water-Ti-Bone-Water\n(config: 3; s_{low}=5 [mm], s_{high}= 10 [mm], voxel = 1 [mm])"
#
set style line 1 lt 1 lw 3 lc rgb 'black'
set style line 2 lt 1 lw 4 lc rgb 'red'
set style line 3 lt 1 dt 4 lw 5 lc rgb 'blue'
set style line 4 lt 1 dt 3 lw 5 lc rgb '#008000'
set style line 5 lt 1 lw 3 lc rgb 'violet'
#
set style line 11 lt 1 pt 7  lw 2 ps 1.2   lc rgb 'black'
set style line 12 lt 1 pt 12  lw 4 ps 1.2    lc rgb 'red'
set style line 13 lt 1 pt 9 lw 2 ps 1.5   lc rgb 'blue'
set style line 14 lt 1 pt 2  lw 3 ps 1    lc rgb '#008000'
set style line 15 lt 1 pt 3  lw 1 ps 1   lc rgb 'violet'
#
set style line 16 lt 1 dt 2  lw 3  lc rgb 'black'
set style line 17 lw 5  lc rgb 'black'
#
set style line 101 lt 1 pt 7  lw 2 ps 0.8    lc rgb 'black'
set style line 102 lt 1 pt 12  lw 3 ps 1.0    lc rgb 'red'
set style line 103 lt 1 pt 9 lw 2 ps 1.1    lc rgb 'blue'
set style line 104 lt 1 pt 2  lw 3 ps 1      lc rgb '#008000'
set style line 105 lt 1 pt 3  lw 1 ps 1      lc rgb 'violet'
set style line 106 lt 1 dt 2  lw 3  lc rgb 'black'
######## Axes settings
set ylabel "Dose [MeV cm2/g]"
set xlabel "Depth [cm]" offset 0,0.5
set grid mxtics xtics lw 1.5 lc rgb 'black', lw 1.5 lc rgb 'gray'
set grid mytics ytics lw 1.5 lc rgb 'black', lw 1.5 lc rgb 'gray'
#
set key Left left bottom box reverse samplen 6 w -7 box opaque
set border back
#
set xrange[:7]
set yrange[0:2.5]
#
#
set arrow 1 from 1,0 to 1,2.5 nohead ls 16
set arrow 2 from 1.5,0 to 1.5,2.5 nohead ls 16
set arrow 3 from 2.5,0 to 2.5,2.5 nohead ls 16
#
#
#############
plot 'DPM_res_Water_Ti_Bone_Water_15MeV_DigiFromPaper_Fig15_DPM.dat'         u 1:2 w l ls 3 title 'DPM (digit. from Fig.15.)',\
     'hist.sim'      u 2:3 w l ls 2 title 'dpm-g4cpp'
unset output

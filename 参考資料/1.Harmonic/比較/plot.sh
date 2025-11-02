

## ＜シミュレーションの結果をグラフにプロット＞ ##
 #!/bin/bash

  T=0.2;

# 平均と理論値をエラーバーも兼ねて gnuplot で画面に出力 #
 gnuplot -e "
   reset;
   set grid;
   set bmargin 4;
   set lmargin 12;
   set key box;
   set xlabel font 'Helvetica,15';
   set ylabel font 'Helvetica,15';
   set tics font 'Helvetica,12';
   set key font 'Helvetica,12';
   set xlabel 'tau [ħωτ]' ;
   set ylabel 'energy [E/ħω]' ;
   plot 0.5/tanh(0.5/$T) title '理論値';
   replot '1th.dat' using 1:2:3 with errorbars title '1次' pt 7 pointsize 1;
   pause -1;
   set ylabel 'specific heat [C/k_B]';
   plot 1.0/( $T*$T*4.0*sinh(0.5/$T)*sinh(0.5/$T) ) title '理論値' ;
   replot '1th.dat' using 1:4:5 with errorbars title '1次' pt 7 pointsize 1;
   replot '6th.dat' using 1:4:5 with errorbars title '高次' pt 7 pointsize 1;
   pause -1;
 "

  exit 0;




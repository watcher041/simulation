

## ＜シミュレーションの結果をグラフにプロット＞ ##
 #!/bin/bash

# 余計なファイルがあるかをチェック #
 if [ ! -e data*.dat ]
   then
       echo -e "\n Error! データファイルが存在しません。\n ";
       exit 1;
   fi
 if [ -e all.dat ]
   then
       rm all.dat;
   fi

# 分割数ごとのデータを１つにまとめる #
 cat data_*.dat > data.dat ;

# 平均値、標準偏差を求める（NR==1の後は1を含めて実行される。） #
 echo -e "\n結果をグラフにプロットします…\n";
 awk '
   NR == 1 {
       before=$1;
       for(i=1;i<NF;i++){
           ave[i]=0.0; var[i]=0.0; all[i]=0;
       }
   }
   {
       if( $1 != before ){

         for(i=1;i<NF;i++){
             ave[i]/=all[i]; var[i]/=all[i]; var[i]-=ave[i]*ave[i];
         }

         printf("%f %f %f %f %f %e %e\n",before,ave[1],sqrt(var[1]),ave[2],sqrt(var[2]),ave[3],sqrt(var[3]));

         for(i=1;i<NF;i++){
             ave[i]=0.0; var[i]=0.0;  all[i]=0;
         }

         before=$1;
       }

       for(i=1;i<NF;i++){
           ave[i]+=$(i+1); var[i]+=($(i+1)*$(i+1)); all[i]++;
       }
   }
   END{
     for(i=1;i<NF;i++){
       ave[i]/=all[i]; var[i]/=all[i]; var[i]-=ave[i]*ave[i];
     }
     printf("%f %f %f %f %f %f \n",before,ave[1],sqrt(var[1]),ave[2],sqrt(var[2]),ave[3],sqrt(var[3]));

   }' data.dat >> all.dat ;

# 平均と理論値をエラーバーも兼ねて gnuplot で画面に出力 #
 gnuplot -e "
   reset;
   set grid;
   set key font 'Helvetica,15';
   set bmargin 4;
   set lmargin 12;
   set xlabel font 'Helvetica,15';
   set ylabel font 'Helvetica,15';
   set xlabel 'T [K]';
   set ylabel 'E/N [K]';
   plot 'all.dat' using 1:2:3 with errorbars title '計算値';
   replot 'energy.dat' w l title '実験値';
   pause -1;
   set ylabel '~C{.8\~}/N [C/k_B]';
   set xtics 1;
   plot 'all.dat' using 1:4:5 with errorbars title '計算値';
   replot 'experiment.dat' w l title '実験値';
   pause -1;
   plot 'all.dat' using 1:6:7 with errorbars title '計算値';
   pause -1;
 "

  exit 0;




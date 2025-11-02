

## ＜シミュレーションの結果をグラフにプロット＞ ##
 #!/bin/bash

# 余計なファイルがあるかをチェック #
 if [ ! -e data*.dat ]
   then
       echo -e "\n Error! データファイルが存在しません。\n ";
       exit 1;
   fi
 if [ -e plot.dat ]
   then
       rm plot.dat;
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
           
         printf("%f %f %f %f %f\n",before,ave[1],sqrt(var[1]),ave[2],sqrt(var[2]));
         
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
     printf("%f %f %f %f %f\n",before,ave[1],sqrt(var[1]),ave[2],sqrt(var[2]));
     
   }' data.dat >> plot.dat ;
 
# 平均と理論値をエラーバーも兼ねて gnuplot で画面に出力 #
 gnuplot -e "
   reset;
   set grid;
   set xrange [0.0:0.51];
   set xlabel 'T';
   set ylabel 'energy\n';
   plot 0.5/tanh(0.5/x) title 'theory';
   replot 'plot.dat' using 1:2:3 with errorbars title 'data';
   pause -1;
   set xlabel 'T';
   set ylabel 'specific';
   plot 1.0/( x*x*4.0*sinh(0.5/x)*sinh(0.5/x) ) title 'theory';
   replot 'plot.dat' using 1:4:5 with errorbars title 'data';
   pause -1;
 "
  
  exit 0;
  



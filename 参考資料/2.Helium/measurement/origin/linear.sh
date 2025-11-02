
 #!/bin/bash

 TEMP=$1;

 # もし、読み取った温度が大きければ、その手前の温度のデータで線形補間する。 #
 awk -v T=${TEMP} '
   { 
      if($1 > T){ 
          y=( ($2-yb)/($1-xb) )*(T-xb)+yb;
          printf("T=%f y=%.7f\n",T,y*0.15);
          exit(0);
      }
      xb=$1; yb=$2;
   }' mass.txt

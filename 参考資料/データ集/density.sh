
# 個-J/K-molの文で割る。（1体あたりの量になる） #

 #!/bin/bash
 
 awk '
 {
     printf("%3d %f %f\n",$1,$2,$3*0.150553525);
 }' change.txt > density.dat


# 個-J/K-molの文で割る。（1体あたりの量になる） #

 #!/bin/bash
 
 awk '
 {
     printf("%f %f\n",$1,$2/8.31445986559);
 }' C.dat > experiment.dat

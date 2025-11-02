
# 個-J/K-molの文で割る。（1体あたりの量になる） #

 #!/bin/bash
 
 awk '
 BEGIN{ all=0.0; }
 {
     if($1>2.10) exit 0;
     all+=$2;
     
 }
 END{printf("%e\n",all);}
 ' experiment.dat

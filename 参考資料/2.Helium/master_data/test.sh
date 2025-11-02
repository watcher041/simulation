
#!/bin/bash

 gnuplot -e "
   reset;
   set grid;
   f(x)=a*x+b;
   a=1;
   b=1;
   fit f(x) 'all.dat' using 1:2 via a,b;
 "
 
 exit 0;

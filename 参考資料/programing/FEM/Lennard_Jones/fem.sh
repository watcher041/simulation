

 #!/bin/bash

 gnuplot -e "
   reset;
   plot 0.5*x title 'classic';
   replot 'fem.dat' using 1:2 w l;
   pause -1;
   plot 0.5 title 'classic';
   replot 'fem.dat' using 1:3 w l;
   pause -1;
 "

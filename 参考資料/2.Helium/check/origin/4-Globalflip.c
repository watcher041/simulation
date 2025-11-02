

 // 「経路全体をひとつのものとしてフリップさせる」 //
  #include "pimc.h"

// グローバルフリップをする関数
  int Globalflip(path *ps)
 {
    /* 変数の定義 */
    int i,n,d,k,count,bn;
    double xi,So,Sn,dS,w,phase,phi,a;
    struct path ps2=(*ps);

    /* 変数の初期化 */
    count=0; bn=0; a=(ps2.ag);

    /* 経路に対してスレッド3つ間隔でフリップ判定を行う */
    for(i=0;i<N;i++){

      /* 経路をフリップする前の重みを計算 */
      (ps2.num)=1; So=0.0; 
      for(n=0;n<3*M;n++){
          So+=( ST(&ps2,&n,&i)+SU(&ps2,&n,&i)+SF(&ps2,&n,&i) );
      }

      /* 経路をフリップさせる */
      for(d=0;d<3;d++){
      
        /* 乱数で定める量の計算 */
        k=(int)round( 3.0*Random() ); phi=Random(); xi=1.0-2.0*Random();
        
        for(n=0;n<3*M;n++){
            phase=2.0*M_PI*(double)k*(( (double)n/(3*M) )+phi);
            (ps2.td[n].pe[i].r[d])=(ps->td[n].pe[i].r[d])+a*xi*cos(phase);
        }
        
      }
      
      /* 力を更新 */
      (ps2.num)=2; 
       for(n=0;n<3*M;n++){
           SF(&ps2,&n,&i);
       }

      /* 経路をフリップした後の重みを計算 */
      (ps2.num)=3; Sn=0.0; 
      for(n=0;n<3*M;n++){
          Sn+=( ST(&ps2,&n,&i)+SU(&ps2,&n,&i)+SF(&ps2,&n,&i) );
      }

      /* フリップさせるかを判定 */
      dS=Sn-So;
      if(dS<0.0){
          w=1.0;
      }else{
          w=exp(-tau*dS);
      }
      if( Random() <= w ){
        for(n=0;n<3*M;n++){
          for(d=0;d<3;d++){
              Periodic(&(ps2.td[n].pe[i].r[d]),&bn);
          } (ps->td[n])=(ps2.td[n]); 
        } count++;
      }else{
        for(n=0;n<3*M;n++){
            (ps2.td[n])=(ps->td[n]);
        }
      }

   }

   return count;
 }



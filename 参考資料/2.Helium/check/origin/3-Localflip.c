

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// ローカルフリップを行う関数
  int Localflip(path *ps)
 {
    /* 変数の定義 */
    int i,n,d,under,count,bn;
    double xi,Sn,So,dS,w,a;
    struct path ps2=(*ps);

    /* 変数の初期化 */
    count=0; bn=0; a=(ps2.al);

    /* 各粒子ごとにフリップ判定を行う */
    for(i=0;i<N;i++){
      for(n=0;n<3*M;n++){

        /* ふたつ下のスレッドを表す */
        under=(n+(3*M-1))%(3*M);

        /* 分子をフリップする前の相対確率を計算 */
        (ps2.num)=1; So=ST(&ps2,&under,&i)+ST(&ps2,&n,&i);
                     So+=SU(&ps2,&n,&i)+SF(&ps2,&n,&i);

        /* 各成分ごとにフリップさせる */
        for(d=0;d<3;d++){
            xi=1.0-2.0*Random();
            (ps2.td[n].pe[i].r[d])=(ps->td[n].pe[i].r[d])+a*xi;
        }
        
        /* 力を更新 */
        (ps2.num)=2; SF(&ps2,&n,&i);

        /* フリップさせた後の重みを計算 */
        (ps2.num)=3; Sn=ST(&ps2,&under,&i)+ST(&ps2,&n,&i);
                     Sn+=SU(&ps2,&n,&i)+SF(&ps2,&n,&i);

        /* 分子をフリップさせるか判定する。 */
        dS=Sn-So;
        if(dS<0.0){
            w=1.0;
        }else{
            w=exp(-tau*dS);
        }
        if( Random() <= w ){     
          for(d=0;d<3;d++){
              Periodic(&(ps2.td[n].pe[i].r[d]),&bn);
          } (ps->td[n])=(ps2.td[n]); count++;
        }else{
           (ps2.td[n])=(ps->td[n]);
       }

     }
   }

   return count;
 }

// 周期境界条件を適用する関数（基本セル内にあるか、距離が近い分子はどれか）
  void Periodic(double *l,int *number)
 {
    static double L_half=0.5*L;
    double delta=L_half*(double)(*number);
    double left=-delta,right=L-delta;

    /* 周期境界条件を適用 */
    while(1){
      if( (*l) < left ){
          (*l) += L;
      }else if( (*l) > right ){
          (*l) -= L;
      }else{
          break;
      }
    }
 } 
 
 

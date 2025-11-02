

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// ローカルフリップを行う関数
  int Localflip(path *ps)
 {
    /* 変数の定義 */
    int i,n,d,under,pn,count,bn;
    double xi,Sn,So,dS,w,a;
    struct path ps2=(*ps);

    /* 変数の初期化 */
    count=0; bn=0; a=(ps2.al);

    /* 各粒子ごとにフリップ判定を行う */
    for(i=0;i<N;i++){
      for(n=0;n<3*M;n++){

        /* 下のスレッドを表す */
        under=(n+(3*M-1))%(3*M);

        /* link先を参照 */
        pn=(ps2.td[under].pe[i].link);

        /* 分子をフリップする前の相対確率を計算 */
        (ps2.num)=1; So=ST(&ps2,&under,&i)+ST(&ps2,&n,&pn);
        So+=SU(&ps2,&n,&pn,&pn)+SF(&ps2,&n,&pn,&pn);

        /* 各成分ごとにフリップさせる */
        for(d=0;d<3;d++){
            xi=1.0-2.0*Random();
            (ps2.td[n].pe[pn].r[d])=(ps->td[n].pe[pn].r[d])+a*xi;
            (ps2.td[n].pe[pn].f[d])=0.0;
        }

        /* フリップさせた後の重みを計算 */
        (ps2.num)=2; Sn=ST(&ps2,&under,&i)+ST(&ps2,&n,&pn);
        Sn+=SU(&ps2,&n,&pn,&pn)+SF(&ps2,&n,&pn,&pn);

        /* 分子をフリップさせるか判定する。 */
        dS=Sn-So;
        if(dS<0.0){
            w=1.0;
        }else{
            w=exp(-tau*dS);
        }

        /* 採択するか判定（構造体のコピーを行う） */
        if( Random() <= w ){
          for(d=0;d<3;d++){
              Periodic(&(ps2.td[n].pe[pn].r[d]),&bn);
          } (ps->td[n])=(ps2.td[n]); (ps->pi[pn])=(ps2.pi[pn]); count++;
        }else{
           (ps2.td[n])=(ps->td[n]); (ps2.pi[pn])=(ps->pi[pn]);
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



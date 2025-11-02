

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// 変位の二乗を求める関数（Cn:現在のスレッドでの番号、Cd:下のスレッド）
  double Tn(path *ps,int *tn,int *pn)
 {
    /* 変数の定義 */
    int d,bn,on,po;
    double answer,dx;

    /* Chin-Actionの係数 */
    constexpr double t1=0.5-T0;
    static double ti[3]={1.0/t1,1.0/t1,0.5/T0};

    /* 変数の初期化 */
    bn=1; on=(*tn+1)%(3*M);
    po=(ps->td[*tn].pe[*pn].link);

    answer=0.0;
    for(d=0;d<3;d++){
        dx=(ps->td[on].pe[po].r[d])-(ps->td[*tn].pe[*pn].r[d]);
        Periodic(&dx,&bn); answer+=dx*dx;
    }
    answer*=ti[(*tn)%3];

    return answer;
 }

// 相互作用を計算する関数（0.全て、1,古い、2.新しい、3,力の更新）
  double Un(path *ps,int *tn,int *pn,int *one)
 {
    /* 変数の定義 */
    int i,j,d,bn;
    double x,dx,total;

    /* Chin-Actionの係数 */
    constexpr double t1=1.0-2.0*T0,v1=1.0/(6.0*t1*t1),v2=1.0-2.0*v1;
    static double vi[3] = {v1,v2,v1};

    /* 値の初期化 */
    total=0.0; bn=1;

    /* ”num”の値で場合分けを行う */
    switch(ps->num){

      /* 全てのポテンシャルを計算 */
      case 0:
        for(i=0;i<N;i++){
          for(j=i+1;j<N;j++){

            x=0.0;
            for(d=0;d<3;d++){
                dx=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[j].r[d]);
                Periodic(&dx,&bn); x+=dx*dx;
            }
            x=sqrt(x); total+=U(&x);

          }
        }
      break;

      /* フリップする前後のポテンシャルを計算 */
      case 1:
      case 2:

        for(i=0;i<N;i++){

          if(i==(*pn)) continue;

          x=0.0;
          for(d=0;d<3;d++){
              dx=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*pn].r[d]);
              Periodic(&dx,&bn); x+=dx*dx;
          }
          x=sqrt(x); total+=U(&x);

        }

        /* 交換を行うときに使用 */
        if( (*pn) == (*one) ) break;

        for(i=0;i<N;i++){

          if( i==(*one) || i==(*pn) ) continue;

          x=0.0;
          for(d=0;d<3;d++){
              dx=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*one].r[d]);
              Periodic(&dx,&bn); x+=dx*dx;
          }
          x=sqrt(x); total+=U(&x);

        }

      break;
    }

    total*=vi[(*tn)%3];

    return total;
 }

// ポテンシャルを計算する関数
  double U(double *dr)
 {
    /* コンパイル時定数の定義 */
    constexpr double v=NUM/RHO;
    constexpr double l=pow(v,1.0/3.0)/RM, rc=0.5*l;
    constexpr int max = (int)round(rc/DIGIT);

    /* 静的変数の定義 */
    static int loop=1;
    static double pot[max];

    /* ポテンシャルの計算 */
    while(loop){

      int i;
      double x,x2,x6,x8,x10,fx,u1,u2,u3,u4,value;

      for(i=1;i<max;i++){

        x=(double)i*DIGIT; x2=x*x; x6=x2*x2*x2; x8=x6*x2; x10=x8*x2;
        u1=A*exp(-ALPHA*x+BETA*x2); u2=C6/x6; u3=C8/x8; u4=C10/x10;

        if(x<D){
            value=( (D/x)-1.0 ); fx=exp(-value*value);
        }else{
            fx=1.0;
        }

        pot[i]=u1-(u2+u3+u4)*fx;

      }

      loop--;

    }

    /* カットオフ判定を行う */
    if( (*dr) < rc ){

      int dist=(int)round( (*dr)/DIGIT );

      return pot[dist];

    }else{

      return 0.0;

    }
 }


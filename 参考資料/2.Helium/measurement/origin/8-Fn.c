

// 「シミュレーションの本体」 //
 #include "pimc.h"

// 相互作用を計算する関数（0.全て、1,古い、2.新しい、3,力の更新）
 double Fn(path *ps,int *tn,int *pn,int *one)
{
   /* 変数の定義 */
   int i,d,bn;
   double eij[3],total,error,x,du,force,value;

   /* Chin-Actionで用いる定数 */
   constexpr double t1 = 1.0-2.0*T0,v1 = 1.0/(6.0*t1*t1);
   static double u0 = ( 1.0-(1.0/t1)+(v1/t1) )/12.0;
   static double ai[3] = {A1,1.0-2.0*A1,A1};

   /* 初期化を行う */
   total=0.0; error=0.0; bn=1;

   /* 力の計算を行う */
   switch(ps->num){

     /* あるスレッドにおける全粒子の力の二乗 */
     case 0:
       for(i=0;i<N;i++){
         for(d=0;d<3;d++){
             value=(ps->td[*tn].pe[i].f[d])*(ps->td[*tn].pe[i].f[d]);
             Kahan(&total,&error,&value);
         }
       }
     break;

     /* フリップ前の力を計算 */
     case 1:

       /* 予め、力の二乗を計算しておく */
       for(i=0;i<N;i++){
         for(d=0;d<3;d++){
             value=(ps->td[*tn].pe[i].f[d])*(ps->td[*tn].pe[i].f[d]);
             Kahan(&total,&error,&value);
         }
       }

       /* フリップ判定の際、変化しない力を残す */
       for(i=0;i<N;i++){

         if(i==(*pn)) continue;

         x=0.0;
         for(d=0;d<3;d++){
             eij[d]=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*pn].r[d]);
             Periodic(&eij[d],&bn); x+=eij[d]*eij[d];
         }
         x=sqrt(x); du=dU(&x);
         for(d=0;d<3;d++){
             force=-du*eij[d]; (ps->td[*tn].pe[i].f[d])-=force;
         }

       }

       /* 交換を行うときに使用 */
       if( (*pn) == (*one) ) break;

       for(i=0;i<N;i++){

         if( i==(*one) || i==(*pn) ) continue;

         x=0.0;
         for(d=0;d<3;d++){
             eij[d]=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*one].r[d]);
             Periodic(&eij[d],&bn); x+=eij[d]*eij[d];
         }
         x=sqrt(x); du=dU(&x);
         for(d=0;d<3;d++){
             force=-du*eij[d]; (ps->td[*tn].pe[i].f[d])-=force;
         }

       }

     break;

     /* フリップした後に変化した力を計算 */
     case 2:

       for(i=0;i<N;i++){

         if(i==(*pn)) continue;

           x=0.0;
           for(d=0;d<3;d++){
               eij[d]=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*pn].r[d]);
               Periodic(&eij[d],&bn); x+=eij[d]*eij[d];
           }
           x=sqrt(x); du=dU(&x);

           /* 交換の判定に用いる（交換を行うときに重複しないために、iの条件も入れる） */
           if( i>(*pn) && (*pn)==(*one) && x < RANGE ){
             if( (ps->pi[*pn].pj[i].nb) == -1 || Random() < 1.0/(double)(3*M) ){
                 (ps->pi[*pn].pj[i].nb)=(*tn);
             }
           }

           for(d=0;d<3;d++){
               force=-du*eij[d];
               (ps->td[*tn].pe[i].f[d])+=force;
               (ps->td[*tn].pe[*pn].f[d])-=force;
           }

         }

         /* 交換を行うときに使用 */
         if( (*pn) != (*one) ){

           for(i=0;i<N;i++){

             if( i==(*pn) || i==(*one) ) continue;

             x=0.0;
             for(d=0;d<3;d++){
                 eij[d]=(ps->td[*tn].pe[i].r[d])-(ps->td[*tn].pe[*one].r[d]);
                 Periodic(&eij[d],&bn); x+=eij[d]*eij[d];
             }
             x=sqrt(x); du=dU(&x);
             for(d=0;d<3;d++){
                 force=-du*eij[d];
                 (ps->td[*tn].pe[i].f[d])+=force;
                 (ps->td[*tn].pe[*one].f[d])-=force;
             }

           }
         }

         /* 力を更新してから、力の二乗を計算する */
         for(i=0;i<N;i++){
           for(d=0;d<3;d++){
               value=(ps->td[*tn].pe[i].f[d])*(ps->td[*tn].pe[i].f[d]);
               Kahan(&total,&error,&value);
           }
         }

       break;

    }

   total*=u0*ai[(*tn)%3];

   return total;
}

// ポテンシャルの微分を計算する関数
  double dU(double *dr)
 {
    /* コンパイル時定数の定義 */
    constexpr double v=(NUM/RHO);
    constexpr double l=pow(v,1.0/3.0)/RM, rc=0.5*l;
    constexpr int max =(int)round(rc/DIGIT);

    /* 静的変数の定義 */
    static int loop=1;
    static double div[max];

    /* ポテンシャルの微分の値を計算 */
    while(loop){

      int i;
      double x,x2,x6,x8,x10,fx,df,u1,u2,u3,u4,value;

      for(i=1;i<max;i++){

          x=(double)i*DIGIT; x2=x*x; x6=x2*x2*x2; x8=x6*x2; x10=x8*x2;
          u1=A*exp(-ALPHA*x+BETA*x2); u2=C6/x6; u3=C8/x8; u4=C10/x10;

        if(x<D){
            value=( (D/x)-1.0 );
            fx=exp(-value*value); df=(2.0*D/x2)*value*fx;
        }else{
            fx=1.0; df=0.0;
        }

        div[i]=(-ALPHA+2.0*BETA*x)*u1+( (6.0*u2+8.0*u3+10.0*u4)*fx/x )-(u2+u3+u4)*df;
        div[i]/=x;

      }

      loop--;

    }

    /* カットオフ判定を行う */
    if( (*dr) < rc ){

      int dist=(int)round( (*dr)/DIGIT );

      return div[dist];

    }else{

      return 0.0;

    }
 }



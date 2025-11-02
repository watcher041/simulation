

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// 期待値を求める関数
  void Sampling(path *ps)
 {
    /* 変数の定義 */
    char dd[20],sd[20];
    int i,n,samp,step,dammy;
    double mca[4],error[4];
    double erg,kin,pot,spec,dE,tp,up,fp,es,ed,pot_plus,value;
    FILE *fd,*fs;

    /* 変数の初期化 */
    dammy=0; (ps->num)=0; pot_plus=Correction(); 

    /* 以前作成したデータがある場合には削除する */
    sprintf(dd,"data_%d.dat",M);
    if ( (fd = fopen(dd, "r")) != NULL ){
        remove(dd); fclose(fd);
    }
    sprintf(sd,"series_%d.dat",M);
    if ( (fs = fopen(sd, "r")) != NULL ){
        remove(sd); fclose(fs);
    }

    /* 実際にサンプリングを行い、期待値を計算 */
    for(samp=1;samp<=SAMP_N;samp++){

      /* 期待値を計算するのに必要な式の初期値（情報落ちを考慮） */
      for(n=0;n<4;n++){
         mca[n]=0.0; error[n]=0.0;
     }

     /* 経路を初期化 */
     Initialize(ps);

     /* 熱化を実施 */
     for(step=1;step<=HEAT_MAX;step++){
         Flipfunc(ps,&step);
      }

      /* 期待値を計算 */
      fs=fopen(sd,"w");
      for(step=1;step<=MCS_MAX;step++){

        /* フリップ判定を行う */
        Flipfunc(ps,&step);

        /* 時系列データを採取 */
        tp=0.0; up=0.0; fp=0.0;
        for(n=0;n<3*M;n++){
          for(i=0;i<N;i++){
              tp+=ST(ps,&n,&i);
          }
          up+=SU(ps,&n,&dammy);
          fp+=SF(ps,&n,&dammy);
        }
        es=-tp+up+3.0*fp;

        /* 時系列データの出力 */
        fprintf(fs,"%d %f\n",step,es);

        /* サンプル値の総和を計算 */
        Kahan(&mca[0],&error[0],&tp);
        Kahan(&mca[1],&error[1],&up);
        Kahan(&mca[2],&error[2],&fp);
        Kahan(&mca[3],&error[3],&es);

      }
      fclose(fs);

      for(n=0;n<4;n++){
         mca[n]/=(double)MCS_MAX;
      }
      kin=(4.5/tau)+( (-mca[0]+mca[2])/(double)(N*M) );
      pot=pot_plus+(mca[1]+2.0*mca[2])/(double)(N*M);

      /* エネルギーのゆらぎを計算 */
      dE=0.0; error[0]=0.0;
      fs=fopen(sd,"r");
      while ( fscanf(fs,"%d %lf\n",&step,&ed) != EOF ){
          value=(ed-mca[3])*(ed-mca[3]);
          Kahan(&dE,&error[0],&value);
      }
      dE/=(double)MCS_MAX;
      fclose(fs);

      /* エネルギー、比熱の計算 */
      erg=kin+pot;
      spec=4.5*(double)M+( (tau*tau*dE-2.0*tau*(mca[0]+3.0*mca[2]))/(double)N );

      /* 計算結果を出力 */
      erg*=EPSILON; kin*=EPSILON; pot*=EPSILON;
      printf(" %2d E=%f [K] K=%f [K] U=%f [K] C=%f \n",samp,erg,kin,pot,spec);
      fd=fopen(dd,"a");
      fprintf(fd,"%f %f %f\n",(tau/EPSILON),erg,spec);
      fclose(fd);

    }
    printf("\n");

    return;
  }

// 総和を計算する関数（Kahanのアルゴリズム）
  void Kahan(double *sum,double *error,double *value)
 {
    double t,y;

    y=(*value)-(*error);      /* 誤差の反映 */
    t=(*sum)+y;               /* 試しに加算 */
    (*error)=( t-(*sum) )-y;  /* 誤差の算出 */
    (*sum)=t;                 /* 加算した値を代入 */
 }
 
// フリップ判定を行う関数
  void Flipfunc(path *ps,int *step)
 {
    Localflip(ps);

    if( (*step)%10==0 ){
        Globalflip(ps);
    }
 }

// ポテンシャルの補正分を計算
  double Correction()
 {
    /* コンパイル時定数の定義 */
    constexpr double v=NUM/RHO;
    constexpr double l=pow(v,1.0/3.0)/RM, rc=0.5*l;

    /* 変数の定義 */
    int i;
    double x,x2,x6,x8,x10,fx,u1,u2,u3,u4,pot,h,value,total;

    /* 変数の初期化 */
    total=0.0; h=(19.0*rc)/(double)MAX; /* rcから10rcまでの積分 */

    /* ポテンシャルの補正項の計算 */
    for(i=0;i<=MAX;i++){

      x=rc+(double)i*h; x2=x*x; x6=x2*x2*x2; x8=x6*x2; x10=x8*x2;
      u1=A*exp(-ALPHA*x+BETA*x2); u2=C6/x6; u3=C8/x8; u4=C10/x10;

      if(x<D){
          value=( (D/x)-1.0 ); fx=exp(-value*value);
      }else{
          fx=1.0;
      }
      pot=x2*(u1-(u2+u3+u4)*fx);

      if(i%MAX==0){
          total+=pot;
      }else{
          total+=2.0*( (i%2)+1.0 )*pot;
      }
      
    }
    total*=2.0*M_PI*( RHO*(RM*RM*RM) )*h/3.0;
    printf(" pot_plus=%e\n",total*EPSILON);
   
    return total;
 }
 
 



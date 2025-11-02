
#include"00_PI.h"

// シミュレーションを実行
void simulation(double *r)
{
    char dd[20];
    int n,samp,step;
    double mca[4];
    double energy,kinetic,potential,specific,dE,tp,up,fp,es,ed;
    FILE *fd,*fs;

    /* 期待値と分布を記録するファイルを作成 */
   sprintf(dd,"data_%d.dat",M);
   fd=fopen(dd,"w");

   /* 実際にサンプリングを行い、期待値を計算 */
   for(samp=1;samp<=N_SAMP;samp++){

     /* 期待値を計算するのに必要な式の初期値 */
     for(n=0;n<4;n++){
         mca[n]=0.0; 
     }

     /* 経路を初期化 */
     reset(r);

     /* 熱化を実施 */
     for(step=1;step<=HEAT_MAX;step++){
         flipfunc(r,step);
     }

     /* 時系列データを作成（比熱の計算にも用いる） */
     sprintf(dd,"series_%d.dat",M);
     fs=fopen(dd,"w");
     
     /* 期待値を計算 */
     for(step=1;step<=MCS_MAX;step++){

       /* フリップ判定を行う */
       flipfunc(r,step);

       /* サンプル値を採取 */
       tp=0.0; up=0.0; fp=0.0;
       for(n=0;n<3*M;n++){
           tp+=Tn(r,n);
           up+=Un(r,n);
           fp+=Fn(r,n);
       }
       es=-tp+up+3.0*fp;

       /* 時系列データの出力 */
       fprintf(fs,"%d %f\n",step,es);

       /* サンプル値の総和を計算 */
        mca[0]+=tp; mca[1]+=up;
        mca[2]+=fp; mca[3]+=es;
     }
     fclose(fs);

     for(n=0;n<4;n++){
         mca[n]/=(double)MCS_MAX;
     }
     kinetic=(1.5/TAU)+( (-mca[0]+mca[2])/(double)M );
     potential=(mca[1]+2.0*mca[2])/(double)M;

     /* エネルギーのゆらぎを計算 */
     dE=0.0;
     sprintf(dd,"series_%d.dat",M);
     fs=fopen(dd,"r");
      while ( fscanf(fs,"%d %lf\n",&step,&ed) != EOF ){
         dE+=(ed-mca[3])*(ed-mca[3]);
     }
     dE/=(double)MCS_MAX;
     fclose(fs);

     /* エネルギー、比熱の計算 */
     energy=kinetic+potential;
     specific=1.5*(double)M+( TAU*TAU*dE-2.0*TAU*(mca[0]+3.0*mca[2]) );

     /* 計算結果を出力 */
     printf(" %2d E=%f K=%f U=%f C=%f\n",samp,energy,kinetic,potential,specific);
     fprintf(fd,"%f %f %f\n",T,energy,specific);
   }
   fclose(fd);
}

// フリップ判定を行う関数
 void flipfunc(double *r,int step)
{
   /* 各スレッドにある分子をフリップさせるか判定する */
   localflip(r);
        
   /* 経路全体をフリップさせるか判定する */
   if(step%10==0){
       globalflip(r);
   }
}
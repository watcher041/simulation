

// 「具体的な計算を行う関数一覧」 //
 #include "PI.h"
 
 void Sampling(double *r)
{
   char dd[20];
   int n,samp,step;
   double mca[3],value[3],error[3];
   double energy,kinetic,potential,specific,dE,tp,up,es,ed;
   FILE *fd,*fs;

    /* 期待値と分布を記録するファイルを作成 */
   sprintf(dd,"data_%d.dat",M);
   fd=fopen(dd,"w");

   /* 実際にサンプリングを行い、期待値を計算 */
   for(samp=1;samp<=N_SAMP;samp++){

     /* 期待値を計算するのに必要な式の初期値（情報落ちを考慮） */
     for(n=0;n<3;n++){
         mca[n]=0.0; error[n]=0.0;
     }

     /* 経路を初期化 */
     Initialize(r);

     /* 熱化を実施 */
     for(step=1;step<=HEAT_MAX;step++){
         Flipfunc(r,step,0);
     }

     /* 時系列データを作成（比熱の計算にも用いる） */
     sprintf(dd,"series_%d.dat",M);
     fs=fopen(dd,"w");
     
     /* 期待値を計算 */
     for(step=1;step<=MCS_MAX;step++){

       /* フリップ判定を行う */
       Flipfunc(r,step,0);

       /* サンプル値を採取 */
       tp=0.0; up=0.0;
       for(n=0;n<M;n++){
           tp+=ST(r,n);
           up+=SU(r,n);
       }
       es=-tp+up;

       /* 時系列データの出力 */
       fprintf(fs,"%d %f\n",step,es);

       /* サンプル値の総和を計算 */
       value[0]=tp; value[1]=up; value[2]=es;
       for(n=0;n<3;n++){
           Kahan(&mca[n],&error[n],value[n]);
       }
     }
     fclose(fs);

     for(n=0;n<3;n++){
         mca[n]/=(double)MCS_MAX;
     }
     kinetic=(0.5/tau)-( mca[0]/(double)M );
     potential=mca[1]/(double)M;

     /* エネルギーのゆらぎを計算 */
     dE=0.0; error[0]=0.0;
     sprintf(dd,"series_%d.dat",M);
     fs=fopen(dd,"r");
      while ( fscanf(fs,"%d %lf\n",&step,&ed) != EOF ){
         value[0]=(ed-mca[2])*(ed-mca[2]);
         Kahan(&dE,&error[0],value[0]);
     }
     dE/=(double)MCS_MAX;
     fclose(fs);

     /* エネルギー、比熱の計算 */
     energy=kinetic+potential;
     specific=0.5*(double)M+( tau*tau*dE-2.0*tau*mca[0]);

     /* 計算結果を出力 */
     printf(" %2d E=%f K=%f U=%f C=%f\n",samp,energy,kinetic,potential,specific);
     fprintf(fd,"%f %f %f\n",tau,energy,specific);
   }
   fclose(fd);

}


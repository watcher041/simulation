
#include "00_header.h"

// 熱化、サンプリングを行う関数
 void sampling(Path *path)
 {
    /* 変数の定義 */
    char dd[20],sd[20];
    int i,n,samp,step;
    double mca[4];
    double erg,kin,pot,spec,dE,tp,up,fp,es,ed,value;
    FILE *fd,*fs;

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

        /* 期待値を計算するのに必要な式の初期値 */
        for(n=0;n<4;n++) mca[n]=0.0;

        /* 経路を初期化 */
        reset(path);

        /* 熱化を実施 */
        for(step=1;step<=HEAT_MAX;step++){
            flipfunc(path,step);
        }

        /* 期待値を計算 */
        fs=fopen(sd,"w");
        for(step=1;step<=MCS_MAX;step++){

            /* フリップ判定を行う */
            flipfunc(path,step);

            /* 時系列データを採取 */
            tp=0.0; up=0.0; fp=0.0;
            for(n=0;n<3*M;n++){
                for(i=0;i<N;i++){
                    tp+=Tn(path,i,n);
                    up+=Un(path,i,n);
                    fp+=Fn(path,i,n);
                }
            }
            for(i=0;i<N;i++){
                if(i != path->pinfo[i].link) printf("%d番目：%d\n",i,path->pinfo[i].link);
            }
            es=-tp+0.5*up+3.0*fp; /* upは二重で計算するため半分にする */

            /* 時系列データの出力 */
            fprintf(fs,"%d %f\n",step,es);

            /* サンプル値の総和を計算 */
            mca[0] += tp; mca[1] += 0.5*up;
            mca[2] += fp; mca[3] += es;

        }
        fclose(fs);

        for(n=0;n<4;n++){
            mca[n]/=(double)MCS_MAX;
        }
        kin=(1.5*(double)RANK/TAU)+( (-mca[0]+mca[2])/(double)(N*M) );
        pot=(mca[1]+2.0*mca[2])/(double)(N*M);

        /* エネルギーのゆらぎを計算 */
        dE=0.0;
        fs=fopen(sd,"r");
        while ( fscanf(fs,"%d %lf\n",&step,&ed) != EOF ){
            value = (ed-mca[3])*(ed-mca[3]);
            dE += value;
        }
        dE/=(double)MCS_MAX;
        fclose(fs);

        /* エネルギー、比熱の計算 */
        erg=kin+pot;
        spec=4.5*(double)M+( (TAU*TAU*dE-2.0*TAU*(mca[0]+3.0*mca[2]))/(double)N );

        /* 計算結果を出力 */
        erg*=EPS; kin*=EPS; pot*=EPS;
        printf(" %2d E=%e [K] K=%e [K] U=%e [K] C=%e \n",samp,erg,kin,pot,spec);
        fd=fopen(dd,"a");
        fprintf(fd,"%f %e %e\n",T*EPS,erg,spec);
        fclose(fd);

        }
        printf("\n");

    return;
  }

// フリップ判定を行う関数
  void flipfunc(Path *path,int step)
 {
    localflip(path);

    if( step%10==0 ){
        globalflip(path);
        
    }
 }
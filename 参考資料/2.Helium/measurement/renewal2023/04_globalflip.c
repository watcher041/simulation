
#include "00_header.h"

// グローバルフリップをする関数
  void globalflip(Path *path)
 {
    /* ループに利用する変数 */
    int i,j,n,d;

    /* フリップに利用する変数 */
    int k;
    double xi,ds,w,phase,phi;
    static double ag = sqrt(ETA*TAU);
    Path path_copy=(*path);

    /* 各粒子ごとにループさせる */
    for(i=0;i<N;i++){

        /* 変更する力を削除しておく */
        for(j=0;j<N;j++){
            if(j==i) continue;
            for(n=0;n<3*M;n++){
                fij(&path_copy,i,j,n,-1);
            }
        }
        for(n=0;n<3*M;n++){
            for(d=0;d<RANK;d++){
                path_copy.pinfo[i].tinfo[n].f[d] = 0.0;
            }
        }

        /* 経路を仮にフリップさせる */
        for(d=0;d<RANK;d++){
        
            /* 乱数で経路全体をフリップさせる */
            k=(int)round( 3.0*Random() ); 
            phi=Random(); xi=1.0-2.0*Random();
            for(n=0;n<3*M;n++){
                phase = 2.0*M_PI*(double)k*(( (double)n/(3*M) )+phi);
                (path_copy.pinfo[i].tinfo[n].r[d]) = (path->pinfo[i].tinfo[n].r[d])+ag*xi*cos(phase);
                periodic(&(path_copy.pinfo[i].tinfo[n].r[d]),2);
            }
            
        }
            
        /* 力を更新 */
        for(j=0;j<N;j++){
            if(j==i) continue;
            for(n=0;n<3*M;n++){
                fij(&path_copy,i,j,n,1);
            }
        }

        /* 経路をフリップする重みを計算 */
        ds=0.0;
        for(n=0;n<3*M;n++){
            ds += Tn(&path_copy,i,n) - Tn(path,i,n);
            ds += Un(&path_copy,i,n) - Un(path,i,n);
            ds += Fn(&path_copy,i,n) - Fn(path,i,n);
        }

        /* フリップさせるかを判定 */
        if(ds<0.0){
            w=1.0;
        }else{
            w=exp(-TAU*ds);
        }
        if( Random() <= w ){
            for(n=0;n<3*M;n++){
                path->pinfo[i].tinfo[n] = path_copy.pinfo[i].tinfo[n];
            } 
        }
        else{
            for(n=0;n<3*M;n++){
                path_copy.pinfo[i].tinfo[n] = path->pinfo[i].tinfo[n];
            }
        }
    }
 }



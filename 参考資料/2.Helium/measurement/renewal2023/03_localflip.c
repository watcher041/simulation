
#include "00_header.h"

// 局所的に経路変更を行う関数
void localflip(Path *path)
{
    /* ループに使用する変数 */
    int i,j,n,d;

    /* スレッドの識別するための変数 */
    int under;

    /* フリップ時の計算に利用するための変数 */
    double ds,w,xi;
    static double al = sqrt(ETA*TAU);

    /* 仮の経路を用意する */
    Path path_copy = (*path);

    /* 各粒子ごとにループする */
    for(i=0;i<N;i++){

        /* スレッドごとにループする */
        for(n=0;n<3*M;n++){
           
            /* 二つ下のスレッドを指定する */
            under = (n+(3*M-1))%(3*M);

            /* フリップ前に変化する力を削除する */
            for(j=0;j<N;j++){
                if(j==i) continue;
                fij(&path_copy,i,j,n,-1);
            }
            for(d=0;d<RANK;d++){
                path_copy.pinfo[i].tinfo[n].f[d] = 0.0;
            }

            /* 経路を仮にフリップさせる */
            for(d=0;d<RANK;d++){
                xi=1.0-2.0*Random();
                (path_copy.pinfo[i].tinfo[n].r[d]) = (path->pinfo[i].tinfo[n].r[d])+al*xi;
                periodic(&(path_copy.pinfo[i].tinfo[n].r[d]),2);
            }

            /* フリップ後の力を計算する */
            for(j=0;j<N;j++){
                if(j==i) continue;
                fij(&path_copy,i,j,n,1);
            }

            /* 分子をフリップする前の相対確率を計算 */
            ds = Tn(&path_copy,i,under)-Tn(path,i,under);
            ds += Tn(&path_copy,i,n)-Tn(path,i,n);
            ds += Un(&path_copy,i,n)-Un(path,i,n);
            ds += Fn(&path_copy,i,n)-Fn(path,i,n);

            /* 経路をフリップさせるか判定する数値を求める */
            if(ds<0.0){
                w = 1.0;
            }
            else{
                w = exp(-TAU*ds);
            }

            /* フリップする場合は、周期境界条件を適用して採用して、
            　　そうでなければ元に戻す */
            if( Random() <= w ){
                (path->pinfo[i].tinfo[n]) = (path_copy.pinfo[i].tinfo[n]);
            }
            else{
                (path_copy.pinfo[i].tinfo[n]) = (path->pinfo[i].tinfo[n]);
            }
        }

    }

    /* 各粒子の交換も行う */
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){

            /* 粒子の交換を行う */
            path_copy.pinfo[i].link = path->pinfo[j].link;
            path_copy.pinfo[j].link = path->pinfo[i].link;

            /* 交換する前後での相対確率を計算 */
            ds  = Tn(&path_copy,i,3*M-1)-Tn(path,i,3*M-1);
            ds += Tn(&path_copy,j,3*M-1)-Tn(path,j,3*M-1);

            /* 経路を交換させるか判定する数値を求める */
            if(ds<0.0){
                w = 1.0;
            }
            else{
                w = exp(-TAU*ds);
            }

            /* 交換する場合は交換、そうでなければ元に戻す */
            if( Random() <= w ){
                (path->pinfo[i].link) = path_copy.pinfo[j].link;
                (path->pinfo[j].link) = path_copy.pinfo[i].link;
            }
            else{
                (path_copy.pinfo[i].link) = path->pinfo[i].link;
                (path_copy.pinfo[j].link) = path->pinfo[j].link;
            }
        }
    }
}
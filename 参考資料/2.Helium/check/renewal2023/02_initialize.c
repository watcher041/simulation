
#include "00_header.h"
#include <time.h>

// 経路を初期化する関数
  void reset(Path *path)
{
    /* ループに利用する変数 */
    int i,j,n,d;

    /* 配置する際に利用する変数 */
    int number;
    double digit;
    const double a = L / (double)NX; /* 原子間距離 */ 

    // 位置、力の初期化
    for(i=0;i<N;i++){

        /* 経路の最端で接続する粒子の番号を初期化 */
        number = i;
        path->pinfo[i].link = i;
        
        /* NX変数で座標を割り振る */
        for(d=0;d<RANK;d++){

            /* 各桁の番号を抽出する */
            digit = number%NX; digit *= a;
            number /= NX;
            for(n=0;n<3*M;n++){    
                path->pinfo[i].tinfo[n].r[d] = digit; /* 粒子の位置を初期化 */
                path->pinfo[i].tinfo[n].f[d] = 0.0;   /* 粒子に働く力を初期化 */
            }
        }
    }

    /* 各粒子に働く力を計算する */
    /* 変化する前の力を消してから変化する分を追加する */
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            for(n=0;n<3*M;n++){
                fij(path,i,j,n,1);
            }
        }
    }
}

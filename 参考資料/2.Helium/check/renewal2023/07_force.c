
#include "00_header.h"

// 粒子間力を計算する関数
  void fij(Path *path,int i,int j,int n,int reversal)
{
    /* ループに使用する変数を定義 */
    int d;
    
    /* 力の計算に使用*/
    double rij,eij[RANK],value,ur;

    /* 今見ている粒子と対になっている粒子に働く力は反対向きに力を加算 */
    rij=0.0;
    for(d=0;d<RANK;d++){

        /* 粒子間の変位を計算 */
        eij[d]=(path->pinfo[i].tinfo[n].r[d])-(path->pinfo[j].tinfo[n].r[d]);

        /* 周期境界条件 */
        periodic(&eij[d],1);
        rij += eij[d]*eij[d];
    }

    /* 力の大きさを計算して代入 */
    /* reversal = -1 であれば、変化前の力を削除
       reversal =  1 であれば、力をそのまま追加する */
    rij = sqrt(rij); ur = dU(rij);
    for(d=0;d<RANK;d++){
        value = -ur*eij[d]/rij;
        value *= (double)reversal;  /* 1:通常通り力を計算、-1:変化前の力のみ削除する */
        (path->pinfo[i].tinfo[n].f[d]) += value;
        (path->pinfo[j].tinfo[n].f[d]) -= value;
    }

}

// ポテンシャルの微分を計算する関数
  double dU(double x)
{ 
    /* 微分の計算に利用 */
    double x2,x6,x8,x10,fx,df,u1,u2,u3,u4,value,du;

    /* 微分の値を計算 */
    x2 = x*x; x6 = x2*x2*x2; x8 = x6*x2; x10 = x8*x2;
    u1 = A*(-ALPHA+2.0*BETA*x)*exp(-ALPHA*x+BETA*x2); 
    u2 = C6/x6; u3 = C8/x8; u4 = C10/x10;
    if(x<D){
        value = ( (D/x)-1.0 );
        fx = exp(-value*value); df = (2.0*D/x2)*value*fx;
    }else{
        fx = 1.0; df = 0.0;
    }
    du = u1;
    du += ( (6.0*u2+8.0*u3+10.0*u4)*fx/x );
    du -= (u2+u3+u4)*df;

    return du;
 }



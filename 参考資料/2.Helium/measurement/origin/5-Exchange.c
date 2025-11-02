

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// 交換を行う関数（１つのスレッドに対して、グローバルフリップを行う）
  int Exchange(path *ps)
 {
    /* 変数の定義 */
    int i,j,d,n,ns,pn,bn0,bn1,count,step;
    int on[THREADS],pi[THREADS],pj[THREADS];
    double So,Sn,dS,w,dt,all,dxi,dxj;
    struct path ps2=(*ps);

    /* τ の分割幅の割合 */
    constexpr double t1=0.5-T0;
    static double tw[3]={t1,t1,2.0*T0};

    /* 変数の初期化 */
    count=0; bn0=0; bn1=1;

    /* 交換判定を行う */
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){

          /* 原子間距離が許容範囲かを判定 */
          ns=(ps2.pi[i].pj[j].nb); if(ns==-1){ continue; }
          (ps->pi[i].pj[j].nb)=-1;

          /* 交換するスレッドの範囲を指定 */
          for(n=0;n<THREADS;n++){
              on[n]=(ns+n)%(3*M);
          }

          /* 距離全体の比を計算（Chin-Actionのため） */
          all=0.0;
          for(n=0;n<THREADS-1;n++){
              all+=tw[on[n]%3];
          }

          /* 粒子番号の識別 */
          pi[0]=i; pj[0]=j;
          for(n=1;n<THREADS;n++){
              pi[n]=(ps2.td[on[n-1]].pe[pi[n-1]].link);
              pj[n]=(ps2.td[on[n-1]].pe[pj[n-1]].link);
          }

          /* 交換する前の重みを計算 */
          So=ST(&ps2,&on[0],&pi[0])+ST(&ps2,&on[0],&pj[0]); (ps2.num)=1;
          for(n=1;n<THREADS-1;n++){
             So+=ST(&ps2,&on[n],&pi[n])+ST(&ps2,&on[n],&pj[n]);
             So+=SU(&ps2,&on[n],&pi[n],&pj[n])+SF(&ps2,&on[n],&pi[n],&pj[n]);
          }

          /* 経路の交換を行う（原子間距離が近いところで交換を行う） */
          pn=pi[0]; pi[0]=pj[0]; pj[0]=pn; 
          (ps2.td[on[0]].pe[pj[0]].link)=(ps->td[on[0]].pe[pi[0]].link); 
          (ps2.td[on[0]].pe[pi[0]].link)=(ps->td[on[0]].pe[pj[0]].link);

          /* 直線の経路を考える */
          for(d=0;d<3;d++){

            /* 交換する経路の幅を計算 */
            dxi=(ps2.td[on[THREADS-1]].pe[pi[THREADS-1]].r[d])-(ps2.td[on[0]].pe[pi[0]].r[d]); Periodic(&dxi,&bn1);
            dxj=(ps2.td[on[THREADS-1]].pe[pj[THREADS-1]].r[d])-(ps2.td[on[0]].pe[pj[0]].r[d]); Periodic(&dxj,&bn1);

            /* 直線の経路を生成 */
            for(n=1;n<THREADS-1;n++){

              /* 経路が直線になるように座標を求めるs */
              dt=tw[on[n-1]%3];
              (ps2.td[on[n]].pe[pi[n]].r[d])=(ps2.td[on[n-1]].pe[pi[n-1]].r[d])+(dt/all)*dxi; Periodic(&(ps2.td[on[n]].pe[pi[n]].r[d]),&bn0);
              (ps2.td[on[n]].pe[pj[n]].r[d])=(ps2.td[on[n-1]].pe[pj[n-1]].r[d])+(dt/all)*dxj; Periodic(&(ps2.td[on[n]].pe[pj[n]].r[d]),&bn0); 

              /* 位置が変化した粒子に関する力を初期化 */
              (ps2.td[on[n]].pe[pi[n]].f[d])=0.0;
              (ps2.td[on[n]].pe[pj[n]].f[d])=0.0;

            }

          }

          /* 設定経路の熱化 */
          for(step=1;step<=ADJUST_MAX;step++){
              Adjustment(&ps2,on,pi); Adjustment(&ps2,on,pj);
          }

          /* 交換後の重みを計算 */
          Sn=ST(&ps2,&on[0],&pj[0])+ST(&ps2,&on[0],&pi[0]); (ps2.num)=2;
          for(n=1;n<THREADS-1;n++){
              Sn+=ST(&ps2,&on[n],&pi[n])+ST(&ps2,&on[n],&pj[n]);
              Sn+=SU(&ps2,&on[n],&pi[n],&pj[n])+SF(&ps2,&on[n],&pi[n],&pj[n]);
          }

          /* Metropolis法で判定 */
          dS=Sn-So;
          if(dS<0.0){
              w=1.0;
          }else{
              w=exp(-tau*dS);
          }
          if( Random() <= w ){
            for(n=0;n<THREADS-1;n++){
                (ps->td[on[n]])=(ps2.td[on[n]]);
            } count++; printf("%d %d %d %e %e\n",i,j,ns,So,Sn);
          }else{
            for(n=0;n<THREADS-1;n++){
                (ps2.td[on[n]])=(ps->td[on[n]]);
            }
          }


      }
    }


    return count;
 }

// 経路ひとつあたりに含まれる粒子数を計算する関数
  double Winding(path *ps)
 {
    int rn[N],i,n,pn,ring;
    double ave;

    for(i=0;i<N;i++){
        rn[i]=1;
    }

    for(i=0;i<N;i++){

      if(rn[i]==0) continue;

      pn=i;
      for(n=0;n<3*M;n++){
            pn=(ps->td[n].pe[pn].link);
      }

      pn=i; 
      do{
        for(n=0;n<3*M;n++){
            pn=(ps->td[n].pe[pn].link);
        } if(pn!=i) rn[pn]--;
      } while(pn!=i);

    }

    ring=0;
    for(i=0;i<N;i++){
        ring+=rn[i];
    }

    ave=(double)N/(double)ring; printf("Winding=%e\n",ave);

    return ave;
 }

// 交換した経路だけ熱化を行う関数
  void Adjustment(path *ps,int *on,int *pn)
 {
    int n,d,under,bn;
    double old[3],dS,w,xi,a;

     a=(ps->al); bn=0;

    for(n=1;n<THREADS-1;n++){

      /* link先を参照 */
      under=n-1;

      /* 分子をフリップする前の相対確率を計算 */
      (ps->num)=1; dS=-( ST(ps,&on[under],&pn[under])+ST(ps,&on[n],&pn[n])+SU(ps,&on[n],&pn[n],&pn[n]) );

       /* 各成分ごとにフリップさせる */
       for(d=0;d<3;d++){
           xi=1.0-2.0*Random(); old[d]=(ps->td[on[n]].pe[pn[n]].r[d]);
           (ps->td[on[n]].pe[pn[n]].r[d])=old[d]+a*xi;
       }

        /* フリップさせた後の重みを計算 */
        (ps->num)=2; dS+=( ST(ps,&on[under],&pn[under])+ST(ps,&on[n],&pn[n])+SU(ps,&on[n],&pn[n],&pn[n]) );

        /* 分子をフリップさせるか判定する。 */
        if(dS<0.0){
            w=1.0;
        }else{
            w=exp(-tau*dS);
        }

        /* 採択するか判定（構造体のコピーを行う） */
        if( Random() <= w ){
          for(d=0;d<3;d++){
              Periodic(&(ps->td[on[n]].pe[pn[n]].r[d]),&bn);
          }
        }else{
          for(d=0;d<3;d++){
              (ps->td[on[n]].pe[pn[n]].r[d])=old[d];
          }
        }

      }
      
      
      
      
 }




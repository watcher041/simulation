

 /* ヘッダーファイルの指定 */
  #include "pimc.h"

// 初期化を行う関数
  void Initialize(path *ps)
 {
    /* 変数の定義 */
    int i,j,n,d,number,nx,bn;
    double eij[3],x,du,force,Lx;

    /* 変数の初期化（Lx:各成分で見たときの原子間距離） */
    bn=1; Lx=L/(double)NL;

    /* 各スレッドごとに位置と力の初期化を行う */
    for(n=0;n<3*M;n++){

      /* 粒子の位置と力の初期値 */
      for(i=0;i<N;i++){

        /* NL進数で初期配置を決定 、力を初期化*/
        number=i;
        for(d=0;d<3;d++){
            nx=(number%NL); number/=NL;
            (ps->td[n].pe[i].r[d]) = (double)nx*Lx;
            (ps->td[n].pe[i].f[d]) = 0.0;
        }

      }

      /* 各スレッドごとに粒子に作用する力を計算する */
      for(i=0;i<N;i++){ 
        for(j=i+1;j<N;j++){

          /* 力の計算（今見ている粒子と対になっている粒子に働く力は反対向きになっている） */
          x=0.0;
          for(d=0;d<3;d++){
              eij[d]=(ps->td[n].pe[i].r[d])-(ps->td[n].pe[j].r[d]);
              Periodic(&eij[d],&bn); x+=eij[d]*eij[d];
          }
          x=sqrt(x); du=dU(&x);
          for(d=0;d<3;d++){
              force=-du*eij[d];
              (ps->td[n].pe[i].f[d])+=force;
              (ps->td[n].pe[j].f[d])-=force;
          }

        }
      }

    }

    return;
 }

// フリップ幅を決める関数
  void Amplitude(path *ps)
 {
    /* 変数の定義 */
    int step,local,global;
    long long int count1,count2,all1,all2;
    double judge;

    /* 変数の初期化 */
     local=0; all1=HEAT_MAX*N*(3*M); (ps->al)=AL;
    global=0; all2=(HEAT_MAX/10)*N;  (ps->ag)=AG;

    /* フリップ割合が5割になるまで、ループする */
    while(local*global==0){

      /* フリップに関わる変数を初期化 */
      local=0; global=0; count1=0; count2=0;

      /* 経路を初期化 */
      Initialize(ps);

      /* 熱化を実施 */
      for(step=1;step<=HEAT_MAX;step++){

        /* 各スレッドにある分子をフリップさせるか判定する */
        count1+=Localflip(ps);

        /* 経路全体をフリップさせるか判定する */
        if(step%10==0){
            count2+=Globalflip(ps);
        }

      }

      /* フリップ割合の計算（ローカルフリップ） */
      judge=(double)count1/(double)all1;
      if(0.46>judge){
          (ps->al)*=0.99;
      }else if(0.54<judge){
          (ps->al)*=1.01;
      }else{
          local=1;
      }
      printf("　al=%e local=%f\n",(ps->al),judge);

      /* フリップ割合の計算（グローバルフリップ） */
      judge=(double)count2/(double)all2;
      if(0.46>judge){
          (ps->ag)*=0.99;
      }else if(0.54<judge){
          (ps->ag)*=1.01;
      }else{
          global=1;
      }
      printf("　ag=%e global=%f\n",(ps->ag),judge);

    }

    /* フリップ幅が決まったら、そのフリップ幅を出力 */
    printf("\n");
    printf("フリップ幅は以下のように決定しました。\n");
    printf("　al=%e\n",(ps->al));
    printf("　ag=%e\n",(ps->ag));
    printf("\n");
    printf("これより、サンプリングに入ります…\n");
    printf("\n");

    return;
 }


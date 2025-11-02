

// 名前空間を介して「Hello world!」を出力するプログラム  //

/* ヘッダの指定 */
  #include <iostream>

/* 名前空間の定義 */
  namespace example {
      
    /* 名前空間（example）の中で関数を定義 */
    void j()
    {
        std::cout << "Hello world!" << std::endl ;
    }
      
  }
  
/* メイン関数 */
  int main()
  {
      /* 名前空間（example）内の関数（j()）を実行 */   
      example::j();
  
      return 0;
  }



// 「Hello world!」を出力するプログラム1 //

/* ヘッダの指定（C言語でのヘッダ”ファイル”と区別） */
  #include <iostream>

/* メイン関数 */
  int main()
  {
      /* ここで出力（ 「std::〜」 はヘッダ内の名前空間”std”にあるものを指す） */
      std::cout << "Hello world!" << std::endl ;
 
      return 0;
  }


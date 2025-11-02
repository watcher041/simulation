

// 「Hello world!」を二回出力するプログラム2 //

/* ヘッダの指定 */
  #include <iostream>
  
/* メイン関数 */
  int main()
  {
      /* 名前空間の指定（最初に書いておくことで、”std::”を省略できる） */
      using std::cout;
      using std::endl;
  
      cout << "Hello world! " << endl << "Hello world!" <<endl;
      
      return 0;
  }
  

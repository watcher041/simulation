
// 環境変数を宣言
include!("00_define.rs");

// シミュレーションの全体像
fn main()
{
    /* シミュレーションの初期値を計算 */
    let beta = 1.0/macro_var!(TEMP); 
    let tau = beta/macro_var!(M) as f64;

    /* 経路を宣言 */ 
    let r: [f64; macro_var!(M)] = [0.0; macro_var!(M)];

    /* 経路を初期化 */
    println!("r={:?}",beta);
    println!("r={:?}",tau);
    println!("r={:?}",r);
}
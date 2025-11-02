
// 環境変数を宣言
include!("define.rs");

// シミュレーションの全体像
fn main()
{
    /* シミュレーションの初期値を計算 */
    let beta = 1.0/macro_var!(TEMP); 
    let tau = beta/macro_var!(M) as f64;

    /* 経路を宣言（同時に一旦初期化が必要） */ 
    let mut r: [f64; macro_var!(M)] = [0.0; macro_var!(M)];

    /* 経路を初期化 */
    path::initialize(&mut r);

    println!("r={:?}",beta);
    println!("r={:?}",tau);
    println!("r={:?}",r);
}
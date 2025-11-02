
// 環境変数を宣言
include!("define.rs");

// シミュレーションの全体像
fn main()
{
    /* シミュレーションの初期値を計算 */
    let beta = 1.0/macrovar!(TEMP); 
    let tau = beta/macrovar!(M) as f64;

    /* 経路を宣言（同時に一旦初期化が必要） */ 
    let mut r: [f64; macrovar!(M)] = [0.0; macrovar!(M)];

    /* 経路を初期化 */
    path::initialize(&mut r);

    /* 乱数の生成方法を初期化 */
    let mut seed = rand::rng();

    /* シミュレーションを実行 */
    path::simulation(&mut seed, &mut r);

    println!("r={:?}",beta);
    println!("r={:?}",tau);
    println!("r={:?}",r);
}
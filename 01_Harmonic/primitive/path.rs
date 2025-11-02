
// 乱数を使うため宣言
use rand::Rng;
use rand::rngs::ThreadRng;

// 経路を初期化する関数
pub fn initialize(r: &mut [f64]) 
{
    for i in 0..r.len() {
        r[i] = 0.0;
    }
}

// 経路に対してシミュレーションを実行する関数
pub fn simulation(seed:&mut ThreadRng, r: &mut [f64])
{
    for i in 0..r.len() {
        r[i] = seed.random();
    }
}
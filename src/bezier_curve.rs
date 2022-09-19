use nalgebra::Vector3;

pub fn factorial(n: usize) -> usize {
    let mut f = 1;
    for i in 1..n + 1 {
        f *= i;
    }
    f
}

pub fn binary_coef(n: usize, k: usize) -> usize {
    factorial(n) / (factorial(n - k) * factorial(k))
}

pub fn bezier_curve(
    p: Vec<Vector3<f64>>,
    t: f64,
) -> Vector3<f64> {
    let n = p.len();
    let mut s = Vector3::zeros();
    for i in 0..n {
        let j = n - 1 - i;
        let c = binary_coef(n - 1, i) as f64;
        s += (1.0 - t).powi(j as i32) * t.powi(i as i32) * p[i] * c;
    }
    s
}
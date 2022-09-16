use nalgebra::{Vector3, DMatrix};



fn factorial(n: usize) -> usize {
    let mut f = 1;
    for i in 1..n + 1 {
        f *= i;
    }
    f
}

fn binary_coef(n: usize, k: usize) -> usize {
    factorial(n) / (factorial(n - k) * factorial(k))
}

pub fn implicit_bezier_curve(trgt_point: Vector3<f64>, ctrl_points: Vec::<Vector3<f64>>) -> Result<f64, &'static str> {
    let n = ctrl_points.len();
    let mut mtx = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let mut line = 0.0;
            let s = i + j + 1;
            for l in 0..s {
                let m = s - l;
                let c = (binary_coef(n, l) * binary_coef(n, m)) as f64;
                line += c * (ctrl_points[l].cross(&ctrl_points[m])).dot(&trgt_point);
            }
            mtx[(i, j)] = line;
        }
    }

    match n {
        2 => Ok(mtx.determinant()),
        3 => { 
            let mut alt = DMatrix::<f64>::zeros(n, n);
            alt[(0, 0)] = mtx.clone().remove_column(2).determinant();
            alt[(0, 1)] = mtx.clone().remove_column(1).determinant();
            alt[(1, 0)] = mtx.clone().remove_column(1).determinant();
            alt[(1, 1)] = mtx.clone().remove_column(0).determinant();
            Ok(alt.determinant())
        },
        _ => Err("Neither n <= 1 nor n >= 4 are implemented.")
    }
}
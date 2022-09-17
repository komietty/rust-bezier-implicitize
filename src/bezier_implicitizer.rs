use nalgebra::{DMatrix, Vector3};

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

pub fn implicit_bezier_curve(ctrl_points: Vec<Vector3<f64>>) -> DMatrix<Vector3<f64>> {
    let n = ctrl_points.len();
    let mut mtx = DMatrix::<Vector3<f64>>::zeros(n - 1, n - 1);
    for i in 0..n - 1 {
        for j in 0..n - 1 {
            let mut line = Vector3::<f64>::zeros();
            let s = i + j + 1;
            for l in 0..n {
                if l > s { continue; }
                let m = s - l;
                if m < n && l < m {
                    let c = (binary_coef(n - 1, l) * binary_coef(n - 1, m)) as f64;
                    line += c * (ctrl_points[l].cross(&ctrl_points[m]));
                    println!(
                        "l: {}, m: {}, c: {}, line: [{}, {}, {}]",
                        l, m, c, line.x, line.y, line.z
                    );
                }
            }
            mtx[(i, j)] = line;
        }
    }
    mtx
}

pub fn calc_deviation(m: &DMatrix<Vector3<f64>>, p: Vector3<f64>) -> Result<f64, &'static str> {
    let n = m.nrows();
    let mut m1 = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            m1[(i, j)] = m[(i, j)].dot(&p);
        }
    }
    match n {
        2 => Ok(m1.determinant()),
        /*
        3 => {
            let mut d = DMatrix::<f64>::zeros(n, n);
            d[(0, 0)] = m1.clone().remove_column(2).determinant();
            d[(0, 1)] = m1.clone().remove_column(1).determinant();
            d[(1, 0)] = m1.clone().remove_column(1).determinant();
            d[(1, 1)] = m1.clone().remove_column(0).determinant();
            Ok(d.determinant())
        },
        */
        _ => Err("Neither n <= 1 nor n >= 4 are implemented."),
    }
}

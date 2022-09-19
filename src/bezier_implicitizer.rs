use nalgebra::{DMatrix, Vector3};
use crate::bezier_curve::{binary_coef};

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
                    /*
                    println!(
                        "i: {}, j: {}, l: {}, m: {}, c: {}, line: [{}, {}, {}]",
                        i, j, l, m, c, line.x, line.y, line.z
                    );
                    */
                }
            }
                    println!(
                        "i: {}, j: {}, line: [{}, {}, {}]",
                        i, j, line.x, line.y, line.z
                    );
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
        3 => {
            /*
            Ok(m1.determinant())
            */
            let mut d = DMatrix::<f64>::zeros(2, 2);
            let mut m3 = DMatrix::<Vector3<f64>>::from_row_slice(2, 3, &[
                Vector3::new(-6.0, 3.0, 0.0),
                Vector3::new(-9.0, 12.0, -9.0),
                Vector3::new(3.0, 9.0, -9.0),

                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(3.0, 3.0, -9.0),
            ]);
            let mut m2 = DMatrix::<f64>::from_row_slice(2, 3, &[
                Vector3::new(-6.0, 3.0, 0.0).dot(&p),
                Vector3::new(-9.0, 12.0, -9.0).dot(&p),
                Vector3::new(3.0, 9.0, -9.0).dot(&p),

                Vector3::new(0.0, 0.0, 0.0).dot(&p),
                Vector3::new(0.0, 3.0, 0.0).dot(&p),
                Vector3::new(3.0, 3.0, -9.0).dot(&p),
            ]);


            d[(0, 0)] = m2.clone().remove_column(2).determinant();
            d[(0, 1)] = m2.clone().remove_column(1).determinant();
            d[(1, 0)] = m2.clone().remove_column(1).determinant();
            d[(1, 1)] = m2.clone().remove_column(0).determinant();
            Ok(d.determinant())
        },
        _ => Err("Neither n <= 1 nor n >= 4 are implemented."),
    }
}

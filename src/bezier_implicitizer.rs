use nalgebra::{DMatrix, Vector3, DVector};
use crate::bezier_curve::{binary_coef};

pub fn implicit_bezier_curve(ctrl_points: &Vec<Vector3<f64>>) -> DMatrix<Vector3<f64>> {
    let n = ctrl_points.len() - 1;
    let mut mtx = DMatrix::<Vector3<f64>>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let mut line = Vector3::<f64>::zeros();
            let s = i + j + 1;
            for l in 0..n {
                if l >= s { continue; }
                let m = s - l;
                if m <= n && l <= i {
                    let c = (binary_coef(n, l) * binary_coef(n, m)) as f64;
                    line += c * (ctrl_points[l].cross(&ctrl_points[m]));
                }
            }
            println!( "i: {}, j: {}, line: [{}, {}, {}]", i, j, line.x, line.y, line.z);
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
            Ok(m1.determinant())
            /*
            let mut d = DMatrix::<f64>::zeros(2, 2);
            let mut m2 = DMatrix::<f64>::from_row_slice(2, 3, &[
                Vector3::new(-6.0, 3.0, 0.0).dot(&p),
                Vector3::new(-9.0, 12.0, -9.0).dot(&p),
                Vector3::new(3.0, 9.0, -9.0).dot(&p),

                //Vector3::new(0.0, 0.0, 0.0).dot(&p),
                //Vector3::new(0.0, 3.0, 0.0).dot(&p),
                //Vector3::new(3.0, 3.0, -9.0).dot(&p),
                Vector3::new(0.0, 3.0, 0.0).dot(&p),
                Vector3::new(3.0, 9.0, -9.0).dot(&p),
                Vector3::new(6.0, 6.0, -18.0).dot(&p),
            ]);


            d[(0, 0)] = m2.clone().remove_column(2).determinant();
            d[(0, 1)] = m2.clone().remove_column(1).determinant();
            d[(1, 0)] = m2.clone().remove_column(1).determinant();
            d[(1, 1)] = m2.clone().remove_column(0).determinant();
            Ok(d.determinant())
            */
        },
        4 => {
            let row0 = m.row(0);
            let row1 = m.row(1);
            let row2 = m.row(2);
            let row3 = m.row(3);

            let m2 = DMatrix::from_row_slice(2, 2, &[
                row0[(0, 0)].x, 
                row1[(0, 0)].x, 
                row0[(0, 0)].y, 
                row1[(0, 0)].y, 
            ]);
            let a = m2.clone().try_inverse().unwrap() * DVector::from_row_slice(&[ row2[(0, 0)].x, row2[(0, 0)].y, ]);
            let b = m2.clone().try_inverse().unwrap() * DVector::from_row_slice(&[ row3[(0, 0)].x, row3[(0, 0)].y, ]);
            let s0 = a[(0, 0)];
            let t0 = a[(1, 0)];
            let s1 = b[(0, 0)];
            let t1 = b[(1, 0)];
            let m3 = DMatrix::from_row_slice(2, 3, &[
                (m[(2, 1)] - m[(0, 1)] * s0 - m[(1, 1)] * t0).dot(&p),
                (m[(2, 2)] - m[(0, 2)] * s0 - m[(1, 2)] * t0).dot(&p),
                (m[(2, 3)] - m[(0, 3)] * s0 - m[(1, 3)] * t0).dot(&p),
                (m[(3, 1)] - m[(0, 1)] * s1 - m[(1, 1)] * t1).dot(&p),
                (m[(3, 2)] - m[(0, 2)] * s1 - m[(1, 2)] * t1).dot(&p),
                (m[(3, 3)] - m[(0, 3)] * s1 - m[(1, 3)] * t1).dot(&p),
            ]);


            let mut d = DMatrix::<f64>::zeros(2, 2);
            d[(0, 0)] = m3.clone().remove_column(2).determinant();
            d[(0, 1)] = m3.clone().remove_column(1).determinant();
            d[(1, 0)] = m3.clone().remove_column(1).determinant();
            d[(1, 1)] = m3.clone().remove_column(0).determinant();
            Ok(d.determinant())
        }
        _ => Err("Neither n <= 1 nor n >= 4 are implemented."),
    }
}

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
            //println!( "i: {}, j: {}, line: [{}, {}, {}]", i, j, line.x, line.y, line.z);
            mtx[(i, j)] = line;
        }
    }
    mtx
}

pub fn calc_deviation(m: &DMatrix<Vector3<f64>>, p: Vector3<f64>) -> Result<f64, &'static str> {
    let n = m.nrows();
    let mut mo = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            mo[(i, j)] = m[(i, j)].dot(&p);
        }
    }
    match n % 2 {
        1 => {
            Ok(mo.determinant())
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
        0 => {
            if n == 2 { 
                return Ok(mo.determinant()); }

            let mut mv = m.clone();
            let itr = m.nrows() / 2 - 1;  

            for k in 0..itr {
                let nr = m.nrows() - (k + 1) * 2;
                let nc = m.ncols() - (k + 1);
                let mi = DMatrix::from_row_slice(2, 2, &[
                    mv.row(0)[(0, 0)].x, mv.row(1)[(0, 0)].x, 
                    mv.row(0)[(0, 0)].y, mv.row(1)[(0, 0)].y, 
                ]).try_inverse().unwrap();

                for i in 2 .. nr + 2 {
                    let r = mv.row(i);
                    let a = &mi * DVector::from_row_slice(&[r[(0, 0)].x, r[(0, 0)].y]);
                    let s = a[(0, 0)];
                    let t = a[(1, 0)];
                    for j in 1 .. nc + 1 {
                        mv[(i, j)] = mv[(i, j)] - mv[(0, j)] * s - mv[(1, j)] * t;
                    }
                }
                mv = mv.remove_fixed_rows::<2>(0).remove_fixed_columns::<1>(0);
            }

            let mut v = vec![];
            let c = n / 2 + 1;
            assert!(mv.nrows() == 2 && mv.ncols() == c);
            for i in 0..2 {
                for j in 0..c {
                    v.push(mv[(i, j)].dot(&p));
                }
            }

            let mr = DMatrix::from_row_slice(2, c, &v);
            let mut d = DMatrix::<f64>::zeros(2, 2);
            d[(0, 0)] = mr.clone().remove_column(2).determinant();
            d[(0, 1)] = mr.clone().remove_column(1).determinant();
            d[(1, 0)] = mr.clone().remove_column(1).determinant();
            d[(1, 1)] = mr.clone().remove_column(0).determinant();
            Ok(d.determinant())
        }
        _ => Err(""),
    }
}

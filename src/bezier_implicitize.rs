use crate::bezier_curve::binary_coef;
use nalgebra::{DMatrix, DVector, Vector3};

pub fn implicitize_bezier(ctrl_points: &Vec<Vector3<f64>>) -> DMatrix<Vector3<f64>> {
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
            mtx[(i, j)] = line;
        }
    }
    mtx
}

pub fn reduce_rank(m: &DMatrix<Vector3<f64>>) -> DMatrix<Vector3<f64>> {
    let mut mm = m.clone();
    match m.nrows() % 2 {
        1 => {
            // Only dim3 case is implemented below. Even though this 
            // seems not working collectly. Might ommit some conditions
            // which are not explicitly mentioned in the thesis.
            let nc = m.ncols() - 1;
            let mi = DMatrix::from_row_slice(2, 2, &[
                    mm.row(0)[(0, 0)].x,
                    mm.row(1)[(0, 0)].x,
                    mm.row(0)[(0, 0)].y,
                    mm.row(1)[(0, 0)].y,
                ]).try_inverse().unwrap();
            let i = 2;
            let r = mm.row(i);
            let a = &mi * DVector::from_row_slice(&[r[(0, 0)].x, r[(0, 0)].y]);
            let s = a[(0, 0)];
            let t = a[(1, 0)];
            for j in 1..nc + 1 {
                mm[(i, j)] = mm[(i, j)] - mm[(0, j)] * s - mm[(1, j)] * t;
            }
            mm[(2, 0)] = Vector3::zeros();
            mm = mm.remove_fixed_rows::<1>(0);
        }
        0 => {
            for k in 0..(m.nrows() / 2 - 1) {
                let nr = m.nrows() - (k + 1) * 2;
                let nc = m.ncols() - (k + 1);
                let mi = DMatrix::from_row_slice(2, 2, &[
                        mm.row(0)[(0, 0)].x,
                        mm.row(1)[(0, 0)].x,
                        mm.row(0)[(0, 0)].y,
                        mm.row(1)[(0, 0)].y,
                    ]).try_inverse().unwrap();

                for i in 2..nr + 2 {
                    let r = mm.row(i);
                    let a = &mi * DVector::from_row_slice(&[r[(0, 0)].x, r[(0, 0)].y]);
                    let s = a[0];
                    let t = a[1];
                    for j in 1..nc + 1 {
                        mm[(i, j)] = mm[(i, j)] - mm[(0, j)] * s - mm[(1, j)] * t;
                    }
                }
                mm = mm.remove_fixed_rows::<2>(0).remove_fixed_columns::<1>(0);
            }
        }
        _ => { }
    }
    mm
}


pub fn calc_deviation(m: &DMatrix<Vector3<f64>>, p: Vector3<f64>) -> Result<f64, &'static str> {
    let (nr, nc) = m.shape();
    if nr == nc {
        let mut m1 = DMatrix::<f64>::zeros(nr, nc);
        for i in 0..nr {
        for j in 0..nc {
            m1[(i, j)] = m[(i, j)].dot(&p);
        }}
        return Ok(m1.determinant());
    } else if nr == 2 && nc == 3 {
        let mut m1 = DMatrix::zeros(nr, nc);
        let mut m2 = DMatrix::<f64>::zeros(nc - 1, nc - 1);
        for i in 0..nr {
        for j in 0..nc {
            m1[(i, j)] = m[(i, j)].dot(&p);
        }}
        m2[(0, 0)] = m1.clone().remove_column(2).determinant();
        m2[(0, 1)] = m1.clone().remove_column(1).determinant();
        m2[(1, 0)] = m1.clone().remove_column(1).determinant();
        m2[(1, 1)] = m1.clone().remove_column(0).determinant();
        Ok(m2.determinant())
    } else {
        Err("Either regular matrix or reduced 2x3 matrix is capable right now.")
    }
}
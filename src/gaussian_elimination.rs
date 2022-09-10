use nalgebra::{DMatrix, DVector};

pub fn gaussian_elimination(mut m: DMatrix<f64>, mut b: DVector<f64>) -> DVector<f64> {
    for i in 0..m.nrows() - 1 {
        for j in (i + 1)..m.nrows() {
            let coef = m[(j, i)] / m[(i, i)];
            m.set_row(j, &(m.row(j) - m.row(i) * coef));
            b[j] -= b[i] * coef;
        }
    }
    for i in (0..m.nrows()).rev() {
        b[i] /= m[(i, i)];
        m.set_row(i, &(m.row(i) / m[(i, i)]));
        for j in 0..i {
            b[j] -= b[i] * m[(j, i)];
            m[(j, i)] = 0.0;
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use crate::{gaussian_elimination::gaussian_elimination};
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn gaussian_elimination_test() {
        let m = DMatrix::from_row_slice(3, 3, &[5.0, -1.0, -1.0, 2.0, 1.0, -3.0, 1.0, 1.0, 1.0]);
        let b = DVector::from_row_slice(&[0.0, -5.0, 6.0]);
        print!("b_self: {} \nb_nalg: {}",
            gaussian_elimination(m.clone(), b.clone()),
            m.try_inverse().unwrap() * b
        );
    }
}

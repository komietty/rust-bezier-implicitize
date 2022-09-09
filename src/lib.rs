
use nalgebra::{DMatrix, DVector};

pub fn gauss_method (mut m: DMatrix<f64>, mut b: DVector<f64>) -> DVector<f64> {
    for i in 0..(m.nrows() - 1) {
        for j in (i + 1)..m.nrows() {
            let coef = m[(j, i)] / m[(i, i)];
            m.set_row(j, &(m.row(j) - m.row(i) * coef));
            b[j] -= b[i] * coef;
        }
    }

    for k in 1 .. m.nrows() + 1 {
        let i = m.nrows() - k;
        b[i] /= m[(i, i)];
        m.set_row(i, &(m.row(i) / m[(i, i)]));
        for j in 0..i {
            b[j] -= b[i] * m[(j, i)];
            m[(j, i)] = 0.0;
        }
    }
    b
}

fn svd() {
}

// bezout's resultant
fn bezout_resultant() {

}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Vector3, DMatrix, DVector};
    use crate::gauss_method;

    #[test]
    fn gauss_method_test() {
        let m = DMatrix::from_row_slice(3, 3, &[
            1.0, 2.0, 4.0,
            3.0, 4.0, 8.0,
            -1.0, 4.0, 2.0,
        ]);
        let b = DVector::from_row_slice(&[7.0, 15.0, 5.0]);
        //let b = gauss_method(m, b);
        //print!("{}", b);
        let m = DMatrix::from_row_slice(3, 3, &[
            5.0, -1.0, -1.0,
            2.0, 1.0, -3.0,
            1.0, 1.0, 1.0,
        ]);
        let b = DVector::from_row_slice(&[0.0, -5.0, 6.0]);
        let b = gauss_method(m, b);
        print!("{}", b);
        //let i = m.try_inverse();
        //print!("{}", i.unwrap() * b);
    }
}

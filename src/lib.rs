use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2};

mod gaussian_elimination;
mod qr_decomposition;

fn svd() {}

// bezout's resultant
pub fn bezout_resultant(x: f64, y: f64, m1: Matrix2x3<f64>, m2: Matrix2x3<f64>) -> f64 {
    let mut m = Matrix2x3::<f64>::zeros();
    m.set_row(0, &(m1.row(0) * x));
    m.set_row(1, &(m1.row(1) * y));
    m = m - m2;
    let mut d = Matrix2::<f64>::zeros();
    d[(0, 0)] = Matrix2::new(m[(0, 0)], m[(0, 1)], m[(1, 0)], m[(1, 1)]).determinant();
    d[(1, 0)] = Matrix2::new(m[(0, 0)], m[(0, 2)], m[(1, 0)], m[(1, 2)]).determinant();
    d[(0, 1)] = Matrix2::new(m[(0, 0)], m[(0, 2)], m[(1, 0)], m[(1, 2)]).determinant();
    d[(1, 1)] = Matrix2::new(m[(0, 1)], m[(0, 2)], m[(1, 1)], m[(1, 2)]).determinant();
    d.determinant()
}

#[cfg(test)]
mod tests {
    use crate::bezout_resultant;
    use nalgebra::{DMatrix, Matrix2x3, Matrix2};

    #[test]
    fn bezout_resultant_test() {
        let x = 3.0;
        let y = 4.5;
        let m1 = Matrix2x3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0);
        let m2 = Matrix2x3::new(2.0, 4.0, 5.0, 3.0, 1.0, 4.0);
        let o = bezout_resultant(x, y, m1, m2);
        println!("b:{}", o);
        let m3 = Matrix2::new(
            5.0  * x - 10.0,
            5.0  * x - y - 7.0,
            5.0  * x - y - 7.0,
            -5.0 * x - 2.0 * y + 11.0);
        println!("b:{}", m3.determinant())
    }
}

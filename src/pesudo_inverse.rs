
#[cfg(test)]
mod tests {
    use nalgebra::{DMatrix, Matrix2x3, Matrix2, SMatrix, Vector3, DVector};

    #[test]
    fn pesudo_inverse_test() {
        let m = DMatrix::from_row_slice(
            10, 3,
            &[
                0.97, 1.86, 0.41, 
                1.23, 2.18, 0.53, 
                0.80, 1.24, 0.62, 
                1.29, 0.98, 0.51, 
                1.10, 1.23, 0.69, 
                0.67, 0.34, 0.54, 
                0.87, 0.26, 0.62, 
                1.10, 0.16, 0.48, 
                1.92, 0.22, 0.71, 
                1.29, 0.12, 0.62, 
            ],
        ); 
        let i = m.clone().pseudo_inverse(1e-8).unwrap();
        let ii = (m.clone().transpose() * m.clone()).try_inverse().unwrap() * m.clone().transpose();
        let a = DVector::from_row_slice(&[
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
            1000.0,
        ]);
        println!("{}", i * a);
        println!("{}", ii);

    }
}
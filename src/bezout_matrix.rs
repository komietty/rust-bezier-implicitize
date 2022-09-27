use nalgebra::{ DMatrix };

pub fn gen_bezout_matrix(m: DMatrix<i128>) -> DMatrix<i128> {
    let n = m.ncols();
    let mut bm = DMatrix::<i128>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let ii = n - 1 - j;
            let min = if i < ii {i} else {ii}; 
            let mut val = 0;
            for k in 0..min + 1 {
                let v1 = j + k;
                let v2 = i - k;
                val += m[(0, v1)] * m[(1, v2)] - m[(0, v2)] * m[(1, v1)];
            }
            bm[(i, j)] = val;
        }
    }
    bm
}

#[test]
fn bezout_matrix_test() {
    let m = DMatrix::from_row_slice(2, 4, &[
        0, -1, 0, 3,
        1, 0, 5, 0 
    ]);
    println!("{}", gen_bezout_matrix(m));
}

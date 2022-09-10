use nalgebra::{DMatrix, DVector};

pub fn qr_decomposition(m: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let nr = m.nrows();
    let nc = m.ncols();
    assert!(nr >= nc);
    let mut q = DMatrix::<f64>::zeros(nr, nr);
    let mut r = DMatrix::<f64>::zeros(nr, nc);
    let mut c = 0;
    for j in 0..nr {
        if j < nc {
            let a = m.column(j);
            let mut q_j = a.clone_owned();
            for i in 0..j {
                let r_ij = a.dot(&q.column(i));
                q_j -= r_ij * q.column(i);
                r[(i, j)] = r_ij;
            }
            let r_jj = q_j.magnitude();
            q_j /= r_jj;
            q.set_column(j, &q_j);
            r[(j, j)] = r_jj;
        } else {
            loop {
                let e = DMatrix::<f64>::identity(nr, nr);
                let a = e.column(c);
                let mut q_c = a.clone_owned();
                for i in 0..nc {
                    q_c -= a.dot(&q.column(i)) * q.column(i);
                }
                if !(q_c == DVector::<f64>::zeros(q_c.nrows())) {
                    q_c /= q_c.magnitude();
                    q.set_column(j, &q_c);
                    break;
                } else {
                    c += 1;
                }
            }
        }
    }
    (q, r)
}

#[cfg(test)]
mod tests {
    use crate::{qr_decomposition::qr_decomposition};
    use nalgebra::{DMatrix};

    #[test]
    fn qr_decomposition_test() {
        let m = DMatrix::from_row_slice(3, 2, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let d1 = m.clone().qr();
        let d2 = qr_decomposition(&m);
        print!("Q: {} \nR: {}", d1.q(), d1.r());
        print!("Q: {} \nR: {}", d2.0, d2.1);
    }
}

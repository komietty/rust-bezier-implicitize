use nalgebra::{ DMatrix, DMatrixSlice, Matrix2, Vector3 };

pub struct Decomposed {
    p: DMatrix<isize>,
    q: DMatrix<isize>,
    b: DMatrix<isize>,
}

pub fn iamin_full(m: DMatrixSlice<isize>) -> (usize, usize) {
    assert!(!m.is_empty(), "The input matrix must not be empty.");
    let mut the_min = unsafe { m.get_unchecked((0, 0)).abs() };
    let mut the_ij = (0, 0);
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            let val = unsafe { m.get_unchecked((i, j)).abs() };
            if val < the_min && val > 0 {
                the_min = val;
                the_ij = (i, j);
            }
        }
    }
    the_ij
}

pub fn is_zeros(m: DMatrixSlice<isize>) -> bool {
    let (nrows, ncols) = m.shape();
    for j in 0..ncols {
        for i in 0..nrows {
            let el = unsafe { m.get_unchecked((i, j)) };
            if *el != 0 {
                return false;
            }
        }
    }
    true
}

pub fn is_mod_zeros(m: DMatrixSlice<isize>, v: isize) -> bool {
    let (nrows, ncols) = m.shape();
    for j in 0..ncols {
        for i in 0..nrows {
            let el = unsafe { m.get_unchecked((i, j)) };
            if *el % v != 0 {
                return false;
            }
        }
    }
    true
}

pub fn modulo_full(m: DMatrixSlice<isize>) -> (usize, usize) {
    assert!(!m.is_empty(), "The input matrix must not be empty.");
    let mut the_val = unsafe { m.get_unchecked((0, 0)) };
    let mut the_ij = (0, 0);
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            let val = unsafe { m.get_unchecked((i, j)) };
            if val % the_val != 0 {
                the_ij = (i, j);
                break;
            }
        }
    }
    the_ij
}

pub fn eij(i: usize, j: usize, n: usize) -> DMatrix<isize> {
    let mut m = DMatrix::identity(n, n);
    m[(i, i)] = 0;
    m[(j, j)] = 0;
    m[(i, j)] = 1;
    m[(j, i)] = 1;
    m
}

pub fn ei(i: usize, n: usize) -> DMatrix<isize> {
    let mut m = DMatrix::identity(n, n);
    m[(i, i)] = -1;
    m
}

pub fn ec(i: usize, j: usize, c: isize, n: usize) -> DMatrix<isize> {
    let mut m = DMatrix::identity(n, n);
    m[(i, j)] = c;
    m
}

fn swap_min(a: &DMatrix<isize>, k: usize) -> Decomposed {
    let s = a.slice((k, k), (a.nrows() - k, a.ncols() - k));
    let m = iamin_full(s);
    let (i, j) = (m.0 + k, m.1 + k);
    let p = eij(k, i, a.nrows());
    let mut q = eij(j, k, a.ncols());
    let mut b = &p * a * &q;
    if b[(k, k)] < 0 {
        let i = ei(k, a.ncols());
        b *= &i;
        q *= &i;
    }
    Decomposed { p: p, q: q, b: b }
}

fn row_null(a: &DMatrix<isize>, k: usize) -> Decomposed {
    let mut b = a.clone();
    let mut p = DMatrix::<isize>::identity(a.nrows(), a.nrows());
    let mut q = DMatrix::<isize>::identity(a.ncols(), a.ncols());
    for i in k + 1..a.nrows() {
        let d = a[(i, k)] / a[(k, k)];
        let e = ec(i, k, -d, a.nrows());
        b = &e * &b;
        p = &e * &p;
    }
    Decomposed { p: p, q: q, b: b }
}

fn col_null(a: &DMatrix<isize>, k: usize) -> Decomposed {
    let mut b = a.clone();
    let mut p = DMatrix::<isize>::identity(a.nrows(), a.nrows());
    let mut q = DMatrix::<isize>::identity(a.ncols(), a.ncols());
    for j in k + 1..a.ncols() {
        let d = a[(k, j)] / a[(k, k)];
        let e = ec(k, j, -d, a.ncols());
        b = &b * &e;
        q = &q * &e;
    }
    Decomposed { p: p, q: q, b: b }
}

fn rem_mod(a: &DMatrix<isize>, k: usize) -> Decomposed {
    let mut p = DMatrix::<isize>::identity(a.nrows(), a.nrows());
    let mut q = DMatrix::<isize>::identity(a.ncols(), a.ncols());
    let idx = modulo_full(a.slice((0, 0), a.shape()));
    let (i, j) = (idx.0 + k, idx.1 + k);
    let d = a[(i, j)] / a[(0, 0)];
    let i1 = ec(i, k, d, a.nrows());
    let i2 = ec(k, j, -1, a.ncols());
    let b = &i1 * a * &i2;
    p *= i1;
    q *= i2;
    Decomposed { p: p, q: q, b: b }
}

pub fn smith_normalize(a: &DMatrix<isize>) -> Decomposed {
    let (nr, nc) = a.shape();
    let mut b = a.clone();
    let mut p = DMatrix::<isize>::identity(nr, nr);
    let mut q = DMatrix::<isize>::identity(nc, nc);
    for k in 0..if nr < nc { nr } else { nc } {
        if is_zeros(b.slice((k, k), (nr - k, nc - k))) { break; }
        loop {
            let d1 = swap_min(&b, k);
            (p, q) = (&d1.p * &p, &q * &d1.q);
            let d2 = row_null(&d1.b, k);
            (p, q) = (&d2.p * &p, &q * &d2.q);

            if is_zeros(d2.b.slice((k + 1, k), (nr - k - 1, 1))) {
                let d3 = col_null(&d2.b, k);
                (p, q) = (&d3.p * &p, &q * &d3.q);
                if is_zeros(d3.b.slice((k, k + 1), (1, nc - k - 1))) {
                    if is_mod_zeros(d3.b.slice((k + 1, k + 1), (nr - k - 1, nc - k - 1)), d3.b[(k, k)]) {
                        b = d3.b;
                        break;
                    } else {
                        let d4 = rem_mod(&d3.b, k);
                        (p, q) = (&d4.p * &p, &q * &d4.q);
                        b = d4.b;
                    }
                } else { b = d3.b; }
            } else { b = d2.b; }
        }
    }
    Decomposed { p: p, q: q, b: b }
}

#[test]
fn smith_swap_test() {
    let m = DMatrix::from_row_slice(3, 4, &[5, 2, 3, 3, 4, 2, -1, 3, 3, 3, 3, 3]);
    let r = swap_min(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
}

#[test]
fn smith_all() {
    let m1 = DMatrix::from_row_slice(2, 3, &[3, 2, 2, 2, 1, 5]);
    let m2 = DMatrix::from_row_slice(3, 4, &[5, 2, 3, 3, 4, 2, -1, 3, 3, 3, 3, 3]);
    let m3 = DMatrix::from_row_slice(3, 4, &[ 2, 0, 0, 0, 0, 2, -6, 8, 0, 3, 6, 4]);
    let r1 = smith_normalize(&m1); 
    let r2 = smith_normalize(&m2); 
    let r3 = smith_normalize(&m3); 
    assert!(&r1.p * m1 * &r1.q == r1.b);
    assert!(&r2.p * m2 * &r2.q == r2.b);
    assert!(&r3.p * m3 * &r3.q == r3.b);
}

#[test]
fn is_zeros_check() {
    let m = DMatrix::from_row_slice(3, 3, &[0, 0, 0, 0, 0, 0, 0, 0, 0 ]);
    assert!(is_zeros(m.slice((0, 0), (0, 0))));
    assert!(is_zeros(m.slice((0, 0), (1, 1))));
    assert!(is_zeros(m.slice((0, 0), (2, 2))));
}

#[test]
fn smith_null_row() {
    let m = DMatrix::from_row_slice(3, 4, &[2, 2, 3, 3, 4, 2, -1, 3, 3, 3, 3, 3]);
    let r = row_null(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
}

#[test]
fn smith_null_col() {
    let m = DMatrix::from_row_slice(3, 4, &[2, 2, -3, 3, 4, 2, -1, 3, 3, 3, 3, 3]);
    let r = col_null(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
}

#[test]
fn smith_rem_modulo() {
    let m = DMatrix::from_row_slice(3, 4, &[2, 0, 0, 0, 0, 2, -6, 8, 0, 3, 6, 4]);
    let r = rem_mod(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
}

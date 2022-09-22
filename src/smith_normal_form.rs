use nalgebra::{DMatrix, Matrix2, SimdPartialOrd, Vector3, Matrix, DMatrixSlice};

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
                if el != &0 { return false; }
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

pub fn eij( i: usize, j: usize, n: usize) -> DMatrix<isize> {
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
    for j in k+1..a.ncols() {
        let d = a[(k, j)] / a[(k, k)];
        let e = ec(k, j, -d, a.ncols());
        b = &b * &e;
        q = &q * &e;
    }
    Decomposed { p: p, q: q, b: b }
}

fn rem_modulo(a: &DMatrix<isize>, k: usize) -> Decomposed {
    let mut p = DMatrix::<isize>::identity(a.nrows(), a.nrows());
    let mut q = DMatrix::<isize>::identity(a.ncols(), a.ncols());
    let idx = modulo_full(a.slice((0, 0), a.shape()));
    let (i, j) = (idx.0 + k, idx.1 + k);
    let d = a[(i, j)] / a[(0, 0)];
    let i1 = ec(i, k,  d, a.nrows());
    let i2 = ec(k, j, -1, a.ncols());
    let b = &i1 * a * &i2;
    p *= i1;
    q *= i2;
    Decomposed { p: p, q: q, b: b }
}

pub fn smith_normalize(a: &DMatrix<isize>) -> Decomposed {
    let (nr, nc) = a.shape();
    let mut a = a.clone();
    let mut p = DMatrix::<isize>::identity(a.nrows(), a.nrows());
    let mut q = DMatrix::<isize>::identity(a.ncols(), a.ncols());
    for k in 0..if nr < nc { nr } else { nc } {
        let slice = a.slice((k, k), (nr - k, nc - k));
        if is_zeros(slice) { break; }
        let mut tmp = 0;
        loop {
            tmp += 1;
            if tmp > 100 { break; }
            
            let d1 = swap_min(&a, k);
            p = &p * &d1.p;
            q = &d1.q * &q;
            println!("b_swaped: {}", d1.b);
            let d2 = row_null(&d1.b, k);
            println!("b_nulled: {}", d2.b);
            println!("slice_2: {}", d2.b.slice((k + 1, k), (nr - k - 1, 1)));
            if is_zeros(d2.b.slice((k + 1, k), (nr - k - 1, 1))) {
                let d3 = col_null(&d2.b, k);
                p = &p * &d3.p;
                q = &d3.q * &q;
                println!("b_col_nulled: {}", d3.b);
                println!("p: {}, dp: {}", p, d3.p);
                println!("q: {}, dq: {}", q, d3.q);
                println!("slice_3: {}", d3.b.slice((k, k + 1), (1, nc - k - 1)));
                if is_zeros(d3.b.slice((k, k + 1), (1, nc - k - 1))) {
                    println!("slice_4: {}", d3.b.slice((k + 1, k + 1), (nr - k - 1, nc - k - 1)));
                    if is_zeros(d3.b.slice((k + 1, k + 1), (nr - k - 1, nc - k - 1))) {
                        return Decomposed { p: p, q: q, b: d3.b };
                    } else {
                        let d4 = rem_modulo(&d3.b, k);
                        p = &p * d4.p;
                        q = &q * d4.q;
                //println!("b_col_nulled: {}", d4.b);
                //println!("p: {}, dp: {}", p, d4.p);
                //println!("q: {}, dq: {}", q, d4.q);
                        panic!("------panipani------");
                    }
                } else {
                    a = d3.b;
                }
            } else {
                a = d1.b;
            }
        }
        /*
        */
    }
    Decomposed { p: p, q: q, b: a }
}

#[test]
fn smith_all() {
    let m = DMatrix::from_row_slice(3, 4, &[
        5, 2, 3, 3,
        4, 2, -1, 3,
        3, 3, 3, 3
        ]
    );
    let m = DMatrix::from_row_slice(2, 3, &[
        3, 2, 2,
        2, 1, 2,
        ]
    );
    let res = smith_normalize(&m);
    //println!("p: {}, q: {}, b:{}", res.p, res.q, res.b);
}

#[test]
fn smith_swap_test() {
    let m = DMatrix::from_row_slice(3, 4, &[
        5, 2, 3, 3,
        4, 2, -1, 3,
        3, 3, 3, 3
        ]
    );
    let r = swap_min(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
}

#[test]
fn smith_null_row() {
    let m = DMatrix::from_row_slice(3, 4, &[
        2, 2, 3, 3,
        4, 2, -1, 3,
        3, 3, 3, 3
        ]
    );
    let r = row_null(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
    println!("p: {}, q: {}, b:{}", r.p, r.q, r.b);
}

#[test]
fn smith_null_col() {
    let m = DMatrix::from_row_slice(3, 4, &[
        2, 2, -3, 3,
        4, 2, -1, 3,
        3, 3, 3, 3
        ]
    );
    let r = col_null(&m, 0);
    assert!(&r.p * m * &r.q == r.b);
    println!("p: {}, q: {}, b:{}", r.p, r.q, r.b);
}

#[test]
fn smith_rem_modulo() {
    let m = DMatrix::from_row_slice(3, 4, &[
        2, 0, 0, 0,
        0, 2, -6, 8,
        0, 3, 6,4 
        ]
    );
    let res = rem_modulo(&m, 0);
    println!("p: {}, q: {}, b:{}", res.p, res.q, res.b);
}
use nalgebra::{DMatrix, Vector3};
use plotters::{prelude::*, style::full_palette::ORANGE_600};
mod bezier_curve;
mod bezier_implicitize;
use crate::{
    bezier_curve::{homogeneous2euclidean,bezier_curve},
    bezier_implicitize::{calc_deviation, implicitize_bezier, reduce_rank},
};

#[test]
fn quartic_test() {
    let circle = |p: &Vector3<f64>, c: RGBColor| Circle::new((p.x as f32, p.y as f32), 2.0, c.filled());
    let r = BitMapBackend::new("quartic.png", (480, 480)).into_drawing_area();
    let mut c = ChartBuilder::on(&r)
        .margin(10)
        .build_cartesian_2d(-1f32..3f32, -2f32..2f32)
        .unwrap();

    r.fill(&WHITE).unwrap();
    c.configure_mesh().draw().unwrap();

    let v: Vec<Vector3<f64>> = vec![
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 1.0, 1.0),
        Vector3::new(1.0, 2.0, 1.0),
        Vector3::new(3.0, 2.0, 1.0),
        Vector3::new(2.0, 0.0, 1.0),
    ];
    let m = reduce_rank(&implicitize_bezier(&v));
    let mut o: Vec<Vector3<f64>> = vec![];

    c.draw_series(LineSeries::new(
        (0..101)
            .map(|t| t as f64 / 100.0)
            .map(|t| bezier_curve(v.clone(), t))
            .map(|hp| homogeneous2euclidean(hp))
            .map(|ep| (ep.x as f32, ep.y as f32)),
        &BLACK,
    )).unwrap();

    for _ in 0..50000 {
        let p = Vector3::<f64>::new(
            rand::random::<f64>() * 3.0,
            rand::random::<f64>() * 4.0 - 2.0,
            1.0,
        );
        if calc_deviation(&m, p).unwrap().abs() < 1.0 { o.push(p); }
    }
    c.draw_series(o.iter().map(|p| circle(p, ORANGE_600))).unwrap();
}

#[test]
fn cubic_test() {
    let circle = |p: &Vector3<f64>, c: RGBColor| Circle::new((p.x as f32, p.y as f32), 2.0, c.filled());
    let r = BitMapBackend::new("cubic.png", (480, 480)).into_drawing_area();
    let mut c = ChartBuilder::on(&r)
        .margin(10)
        .build_cartesian_2d(-1f32..5f32, -3f32..3f32)
        .unwrap();

    r.fill(&WHITE).unwrap();
    c.configure_mesh().draw().unwrap();

    let v: Vec<Vector3<f64>> = vec![
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 1.0, 1.0),
        Vector3::new(1.0, 2.0, 1.0),
        Vector3::new(3.0, 0.0, 1.0),
    ];
    let mut o: Vec<Vector3<f64>> = vec![];
    //let m = reduce_rank(&implicitize_bezier(&v));
    let m = &implicitize_bezier(&v);
    c.draw_series(LineSeries::new(
            (0..101)
                .map(|t| t as f64 / 100.0)
                .map(|t| bezier_curve(v.clone(), t))
                .map(|hp| homogeneous2euclidean(hp))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &BLACK,
        ))
        .unwrap();
    
    for _ in 0..50000 {
        let p = Vector3::<f64>::new(
            rand::random::<f64>() * 6.0,
            rand::random::<f64>() * 6.0 - 3.0,
            1.0,
        );
        if calc_deviation(&m, p).unwrap().abs() < 1.0 { o.push(p); }
    }
    c.draw_series(o.iter().map(|p| circle(p, ORANGE_600))).unwrap();
}

#[test]
fn bezout_matrix_test() {
    fn gen_bezout_matrix(m: DMatrix<f64>) -> DMatrix<f64> {
        let n = m.ncols();
        let mut bm = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let ii = n - 1 - j;
                let min = if i < ii {i} else {ii}; 
                let mut val = 0.0;
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
    let m = DMatrix::from_row_slice(2, 4, &[ 0.0, -1.0, 0.0, 3.0, 1.0, 0.0, 5.0, 0.0 ]);
    println!("bezout matrix: {}, ", gen_bezout_matrix(m));
}
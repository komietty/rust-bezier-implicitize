use nalgebra::Vector3;
use plotters::prelude::*;

use crate::{
    bezier_curve::bezier_curve,
    bezier_implicitizer::{implicitize_bezier, calc_deviation, reduce_rank},
    homogeneous2euclidean,
};

#[test]
fn quartic_test() {
    let r = BitMapBackend::new("quartic.png", (480, 480)).into_drawing_area();
    let mut c = ChartBuilder::on(&r)
        .margin(10)
        .build_cartesian_2d(-1f32..3f32, -2f32..2f32)
        .unwrap();

    r.fill(&RGBColor(41, 45, 62)).unwrap();
    c.configure_mesh().draw().unwrap();

    let v: Vec<Vector3<f64>> = vec![
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 1.0, 1.0),
        Vector3::new(1.0, 2.0, 1.0),
        Vector3::new(3.0, 2.0, 1.0),
        Vector3::new(2.0, 0.0, 1.0),
        //Vector3::new(0.0, 0.0, 1.0),
        //Vector3::new(0.0, 1.0, 1.0),
        //Vector3::new(0.5, 1.5, 1.0),
        //Vector3::new(1.0, 2.0, 1.0),
        //Vector3::new(3.0, 2.0, 1.0),
        //Vector3::new(3.0, 1.0, 1.0),
        //Vector3::new(2.0, 0.0, 1.0),
    ];
    let m = reduce_rank(&implicitize_bezier(&v));
    let mut o: Vec<Vector3<f64>> = vec![];

    c.draw_series(
        v.iter()
            .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, WHITE.filled())),
    )
    .unwrap();

    c.draw_series(LineSeries::new(
        (0..101)
            .map(|t| t as f64 / 100.0)
            .map(|t| bezier_curve(v.clone(), t))
            .map(|hp| homogeneous2euclidean(hp))
            .map(|ep| (ep.x as f32, ep.y as f32)),
        &WHITE,
    ))
    .unwrap();

    for _ in 0..50000 {
        let p = Vector3::<f64>::new(
            rand::random::<f64>() * 3.0,
            rand::random::<f64>() * 4.0 - 2.0,
            1.0,
        );
        if calc_deviation(&m, p).unwrap().abs() < 1.0 {
            o.push(p);
        }
    }

    c.draw_series(
        o.iter()
            .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, GREEN.filled())),
    )
    .unwrap();
}

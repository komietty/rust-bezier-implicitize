use nalgebra::Vector3;
use plotters::prelude::*;

use crate::{
    bezier_curve::bezier_curve,
    bezier_implicitizer::{calc_deviation, implicit_bezier_curve},
    homogeneous2euclidean,
};

#[test]
fn quartic_test() {
    let r = BitMapBackend::new("quartic.png", (480, 480)).into_drawing_area();
    let mut c = ChartBuilder::on(&r)
        .margin(10)
        .build_cartesian_2d(-1f32..4f32, -1f32..4f32)
        .unwrap();

    r.fill(&RGBColor(41, 45, 62)).unwrap();
    c.configure_mesh().draw().unwrap();

    let points: Vec<Vector3<f64>> = vec![
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 1.0, 1.0),
        Vector3::new(1.0, 2.0, 1.0),
        Vector3::new(3.0, 2.0, 1.0),
        Vector3::new(2.0, 0.0, 1.0),
    ];

    c.draw_series(
        points
            .iter()
            .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, WHITE.filled())),
    )
    .unwrap();

    c.draw_series(LineSeries::new(
        (0..101)
            .map(|t| t as f64 / 100.0)
            .map(|t| bezier_curve(points.clone(), t))
            .map(|hp| homogeneous2euclidean(hp))
            .map(|ep| (ep.x as f32, ep.y as f32)),
        &WHITE,
    ))
    .unwrap();

    let m = implicit_bezier_curve(points);
}

use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2, Vector2, Vector3, Vector4, Matrix4};
use plotters::prelude::*;
use rand;
use crate::bezier_implicitizer::{calc_deviation, implicitize_bezier, reduce_rank};
use crate::bezier_curve::{bezier_curve, homogeneous2euclidean, };

#[test]
fn cubic_test() {
    let bg = RGBColor(41, 45, 62);
    let root = BitMapBackend::new("cubic.png", (480, 480)).into_drawing_area();
    root.fill(&bg).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .build_cartesian_2d(-1f32..5f32, -3f32..3f32)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    let p0 = Vector3::new(0.0, 0.0, 1.0);
    let p1 = Vector3::new(0.0, 1.0, 1.0);
    let p2 = Vector3::new(1.0, 2.0, 1.0);
    let p3 = Vector3::new(3.0, 0.0, 1.0);
    let mut points: Vec<Vector3<f64>> = vec![];
    let m1 = reduce_rank(&implicitize_bezier(&vec![p0, p1, p2, p3]));

    chart
        .draw_series(
            vec![
                Circle::new((p0.x as f32, p0.y as f32), 2.0, WHITE.filled()),
                Circle::new((p1.x as f32, p1.y as f32), 2.0, WHITE.filled()),
                Circle::new((p2.x as f32, p2.y as f32), 2.0, WHITE.filled()),
                Circle::new((p3.x as f32, p3.y as f32), 2.0, WHITE.filled()),
            ]
            .into_iter(),
        )
        .unwrap();
    chart
        .draw_series(LineSeries::new(
            (0..101)
                .map(|t| t as f64 / 100.0)
                .map(|t| bezier_curve(vec![p0, p1, p2, p3], t))
                .map(|hp| homogeneous2euclidean(hp))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &WHITE,
        ))
        .unwrap();
    
    for _ in 0..100000 {
        let p = Vector3::<f64>::new(
            rand::random::<f64>() * 6.0,
            rand::random::<f64>() * 6.0-3.0,
            1.0,
        );
        if calc_deviation(&m1, p).unwrap().abs() < 0.5 { points.push(p); }
    }
    chart
        .draw_series(
            points
                .iter()
                .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, GREEN.filled())),
        )
        .unwrap();
}

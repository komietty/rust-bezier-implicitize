use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2, Vector2, Vector3};
use plotters::prelude::*;

mod gaussian_elimination;
mod pesudo_inverse;
mod qr_decomposition;

pub fn homogeneous2euclidean(h: Vector3<f64>) -> Vector2<f64> {
    assert!(h.z != 0.0);
    //Vector2::new(h.x / h.z, h.y / h.z)
    Vector2::new(h.x / h.z, h.y / h.z)
}

pub fn quadratic_bezier_curve(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
) -> Vector3<f64> {
    p0 * (1.0 - t) * (1.0 - t) + p1 * 2.0 * (1.0 - t) * t + p2 * t * t
}

pub fn quadratic_bezier_curve_by_line_in_pencil(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
    k: f64
) -> Vector3<f64> {
    if      k == 0.0 { p0.cross(&(2.0 * p1 * (1.0 - t) + p2 * t)) }
    else if k == 1.0 { p2.cross(&(2.0 * p1 * t + p0 * (1.0 - t))) }
    else { unimplemented!() }
}

pub fn quadratic_bezier_curve_by_lines_of_pencils(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
) -> Vector3<f64> {
    let l0 = p0.cross(&(2.0 * p1 * (1.0 - t) + p2 * t));
    let l1 = p2.cross(&(2.0 * p1 * t + p0 * (1.0 - t)));
    l0.cross(&l1)
}

#[test]
fn chart_context() {
    let bg = RGBColor(41, 45, 62);
    let root = BitMapBackend::new("chart.png", (480, 480)).into_drawing_area();
    root.fill(&bg).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .build_cartesian_2d(-1f32..1f32, -1f32..1f32)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    let p0 = Vector3::new(-0.5, -0.5, 1.0);
    let p1 = Vector3::new(0.0, 0.0, 1.0);
    let p2 = Vector3::new(0.5, -0.2, 1.0);

    chart
        .draw_series(vec![
                Circle::new((p0.x as f32, p0.y as f32), 2.0, WHITE.filled()),
                Circle::new((p1.x as f32, p1.y as f32), 2.0, WHITE.filled()),
                Circle::new((p2.x as f32, p2.y as f32), 2.0, WHITE.filled())
        ].into_iter()).unwrap();

    for t in 0..2 {
        let l0 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, t as f64, 0.0);
        let l1 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, t as f64, 1.0);
        let it = (-100..100) .map(|x| x as f64 * 0.01);
        chart
            .draw_series(LineSeries::new(
                it.clone().map(|x| (x as f32, (-l0.x / l0.y *  x - l0.z / l0.y) as f32)),
                &GREEN,
            ))
            .unwrap();
        chart
            .draw_series(LineSeries::new(
                it.clone().map(|x| (x as f32, (-l1.x / l1.y *  x - l1.z / l1.y + 1e-2) as f32)),
                &BLUE,
            ))
            .unwrap();
    }

    chart
        .draw_series(LineSeries::new(
            (0..100)
                .map(|t| t as f64 / 100.0)
                .map(|t| quadratic_bezier_curve_by_lines_of_pencils(p0, p1, p2, t))
                .map(|hp| homogeneous2euclidean(hp))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &WHITE,
        ))
        .unwrap();

    /*
    chart
        .draw_series(
            (0..100)
                .map(|t| t as f32 * 0.01)
                .map(|t| {
                    let x = t;
                    let y = -t;
                    (x, y)
                })
                .map(|(x, y)| Circle::new((x, y), 1.0, WHITE.filled())),
        )
        .unwrap();
    */
}
use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2, Vector2, Vector3, Vector4, Matrix4};
use plotters::prelude::*;
use rand;
use crate::bezier_implicitizer::{calc_deviation, implicit_bezier_curve};
mod bezier_implicitizer;
mod bezier_curve;
mod bezier_test_quartic;
mod smith_normal_form;

pub fn homogeneous2euclidean(h: Vector3<f64>) -> Vector2<f64> {
    assert!(h.z != 0.0);
    Vector2::new(h.x / h.z, h.y / h.z)
}

pub fn quadratic_bezier_curve_by_line_in_pencil(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
    k: f64,
) -> Vector3<f64> {
    if k == 0.0 {
        p0.cross(&(2.0 * p1 * (1.0 - t) + p2 * t))
    } else if k == 1.0 {
        p2.cross(&(2.0 * p1 * t + p0 * (1.0 - t)))
    } else {
        unimplemented!()
    }
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
    let m1 = implicit_bezier_curve(&vec![p0, p1, p2, p3]);

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
                .map(|t| bezier_curve::bezier_curve(vec![p0, p1, p2, p3], t))
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

#[test]
fn quadratic_test() {
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

    let p3 = Vector3::new(0.5, 0.5, 1.0);
    let p4 = Vector3::new(0.0, -1.0, 1.0);
    let p5 = Vector3::new(-0.5, 0.2, 1.0);

    chart
        .draw_series(
            vec![
                Circle::new((p0.x as f32, p0.y as f32), 2.0, WHITE.filled()),
                Circle::new((p1.x as f32, p1.y as f32), 2.0, WHITE.filled()),
                Circle::new((p2.x as f32, p2.y as f32), 2.0, WHITE.filled()),
            ]
            .into_iter(),
        )
        .unwrap();

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

    chart
        .draw_series(LineSeries::new(
            (0..100)
                .map(|t| t as f64 / 100.0)
                .map(|t| quadratic_bezier_curve_by_lines_of_pencils(p3, p4, p5, t))
                .map(|hp| homogeneous2euclidean(hp))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &WHITE,
        ))
        .unwrap();

    for t in 0..2 {
        let l0 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, t as f64, 0.0);
        let l1 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, t as f64, 1.0);
        let it = (-100..100).map(|x| x as f64 * 0.01);
        chart
            .draw_series(LineSeries::new(
                it.clone()
                    .map(|x| (x as f32, (-l0.x / l0.y * x - l0.z / l0.y) as f32)),
                &GREEN,
            ))
            .unwrap();
        chart
            .draw_series(LineSeries::new(
                it.clone()
                    .map(|x| (x as f32, (-l1.x / l1.y * x - l1.z / l1.y + 1e-2) as f32)),
                &BLUE,
            ))
            .unwrap();
    }

    let l00 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, 0.0, 0.0);
    let l01 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, 1.0, 0.0);
    let l10 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, 0.0, 1.0);
    let l11 = quadratic_bezier_curve_by_line_in_pencil(p0, p1, p2, 1.0, 1.0);
    let s00 = quadratic_bezier_curve_by_line_in_pencil(p3, p4, p5, 0.0, 0.0);
    let s01 = quadratic_bezier_curve_by_line_in_pencil(p3, p4, p5, 1.0, 0.0);
    let s10 = quadratic_bezier_curve_by_line_in_pencil(p3, p4, p5, 0.0, 1.0);
    let s11 = quadratic_bezier_curve_by_line_in_pencil(p3, p4, p5, 1.0, 1.0);

    print!("lines: {}, {}, {}, {}", l00, l01, l10, l11);
    let detA = |x: Vector3<f64>| l00.dot(&x) * l11.dot(&x) - l01.dot(&x) * l10.dot(&x);
    let detB = |x: Vector3<f64>| s00.dot(&x) * s11.dot(&x) - s01.dot(&x) * s10.dot(&x);
    for i in -100..100 {
        let t = i as f64 * 0.1;
        assert!(detA(bezier_curve::bezier_curve(vec![p0, p1, p2], t)).abs() < 1e-10);
    }

    let m1 = implicit_bezier_curve(&vec![p0, p1, p2]);
    let m2 = implicit_bezier_curve(&vec![p3, p4, p5]);

    {
        let mut points: Vec<Vector3<f64>> = vec![];
        let mut points_: Vec<Vector3<f64>> = vec![];
        for _ in 0..10000 {
            let p = Vector3::<f64>::new(
                rand::random::<f64>() * 2.0 - 1.0,
                rand::random::<f64>() * 2.0 - 1.0,
                1.0,
            );
            if detA(p).abs() < 0.05 && detB(p).abs() < 0.05 {
                points.push(p);
            }
            if 
            calc_deviation(&m1, p).unwrap().abs() < 0.02 
            //&& 
            //calc_deviation(&m2, p).unwrap().abs() < 0.02
             {
                points_.push(p);
            }
        }
        chart
            .draw_series(
                points
                    .iter()
                    .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, GREEN.filled())),
            )
            .unwrap();
        chart
            .draw_series(
                points_
                    .iter()
                    .map(|p| Circle::new((p.x as f32, p.y as f32), 2.0, RED.filled())),
            )
            .unwrap();
    }
}

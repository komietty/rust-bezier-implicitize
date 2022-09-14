use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2, Vector2, Vector3};
use plotters::prelude::*;

mod gaussian_elimination;
mod pesudo_inverse;
mod qr_decomposition;

pub fn homogeneous2euclidean(h: Vector3<f64>) -> Vector2<f64> {
    assert!(h.z != 0.0);
    Vector2::new(h.x / h.z, h.y / h.z)
}

// bezout's resultant
pub fn bezout_resultant(x: f64, y: f64, m1: Matrix2x3<f64>, m2: Matrix2x3<f64>) -> f64 {
    let mut m = Matrix2x3::<f64>::zeros();
    m.set_row(0, &(m1.row(0) * x));
    m.set_row(1, &(m1.row(1) * y));
    m = m - m2;
    let mut d = Matrix2::<f64>::zeros();
    d[(0, 0)] = Matrix2::new(m[(0, 0)], m[(0, 1)], m[(1, 0)], m[(1, 1)]).determinant();
    d[(1, 0)] = Matrix2::new(m[(0, 0)], m[(0, 2)], m[(1, 0)], m[(1, 2)]).determinant();
    d[(0, 1)] = Matrix2::new(m[(0, 0)], m[(0, 2)], m[(1, 0)], m[(1, 2)]).determinant();
    d[(1, 1)] = Matrix2::new(m[(0, 1)], m[(0, 2)], m[(1, 1)], m[(1, 2)]).determinant();
    d.determinant()
}

#[cfg(test)]
mod tests {
    use crate::bezout_resultant;
    use nalgebra::{DMatrix, Matrix2, Matrix2x3};

    #[test]
    fn bezout_resultant_test() {
        let x = 3.0;
        let y = 4.5;
        let m1 = Matrix2x3::new(1.0, 2.0, 3.0, 1.0, 2.0, 3.0);
        let m2 = Matrix2x3::new(2.0, 4.0, 5.0, 3.0, 1.0, 4.0);
        let o = bezout_resultant(x, y, m1, m2);
        println!("b:{}", o);
        let m3 = Matrix2::new(
            5.0 * x - 10.0,
            5.0 * x - y - 7.0,
            5.0 * x - y - 7.0,
            -5.0 * x - 2.0 * y + 11.0,
        );
        println!("b:{}", m3.determinant())
    }
}

pub fn quadratic_bezier_curve(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
) -> Vector3<f64> {
    p0 * (1.0 - t) * (1.0 - t) + p1 * 2.0 * (1.0 - t) * t + p2 * t * t
}

pub fn quadratic_bezier_curve_by_lines_of_pencil_0(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
) -> Vector3<f64> {
    let l0 = p0.cross(&(2.0 * p1 * (1.0 - t) + p2 * t));
    l0
}

pub fn quadratic_bezier_curve_by_lines_of_pencils(
    p0: Vector3<f64>,
    p1: Vector3<f64>,
    p2: Vector3<f64>,
    t: f64,
) -> Vector3<f64> {
    let l0 = p0.cross(&(2.0 * p1 * (1.0 - t) + p2 * t));
    let l1 = p2.cross(&(2.0 * p1 * t + p0 * (1.0 - t)));
    l0 * t + l1 * (1.0 - t)
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

    let p0 = Vector3::new(0.2, -0.8, 1.0);
    let p1 = Vector3::new(0.0, 0.0, 1.0);
    let p2 = Vector3::new(0.8, 0.6, 1.0);

    /*
    for t in 0..1 {
        let line = quadratic_bezier_curve_by_lines_of_pencil_0(p0, p1, p2, t as f64);
        let coef = homogeneous2euclidean(line);
        chart
            .draw_series(LineSeries::new(
                (-1..2)
                .map(|t| t as f64)
                .map(|x| (x as f32, (coef.x * x + coef.y) as f32)),
                &GREEN,
            ))
            .unwrap();
    }
    */

    chart
        .draw_series(LineSeries::new(
            (0..100)
                .map(|t| t as f64 / 100.0)
                .map(|t| quadratic_bezier_curve(p0, p1, p2, t))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &BLUE,
        ))
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..100)
                .map(|t| t as f64 / 100.0)
                .map(|t| quadratic_bezier_curve_by_lines_of_pencils(p0, p1, p2, t))
                //.map(|hp| homogeneous2euclidean(hp))
                .map(|ep| (ep.x as f32, ep.y as f32)),
            &RED,
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
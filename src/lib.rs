use crate::{
    bezier_curve::homogeneous2euclidean,
    bezier_implicitizer::{calc_deviation, implicitize_bezier},
};
use nalgebra::{DMatrix, Matrix2, Matrix2x3, Matrix3x2, Matrix4, Vector2, Vector3, Vector4};
use plotters::prelude::*;
mod bezier_curve;
mod bezier_implicitizer;
mod bezout_matrix;
mod test_cubic;
mod test_quartic;
mod test_quadratic;

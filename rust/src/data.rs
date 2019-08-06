// use num::traits::NumOps;
use crate::common;
// use std::ops;

// pub fn squared_difference<T: ops::Mul<Output = T> + ops::Sub<Output = T> + Copy>(a: T, b: T) -> T {
//     (a - b) * (a - b)

pub fn absolute_difference_u8_to_f32(a: u8, b: u8) -> f32 {
    (f32::from(a) - f32::from(b)).abs()
}

pub fn truncated_linear_f32_fh(lambda: f32, tau: f32, a: u8, b: u8) -> f32 {
    lambda * f32::min(tau, absolute_difference_u8_to_f32(a, b))
}

// pub fn truncated_linear_FH<T>(lambda : T, tau : T, a : &u8, b : &u8) -> T where T: NumOps + PartialOrd {
//     // lambda * T::min(tau, (num::FromPrimitive::from_u8(*a).unwrap() - num::FromPrimitive::from_u8(*b).unwrap()).abs());
//     let ab_abs : T = (num::FromPrimitive::from_u8(*a).unwrap() - num::FromPrimitive::from_u8(*b).unwrap()).abs();
//     let this_min = if tau < ab_abs {
//         tau
//     } else {
//         ab_abs
//     };
//     lambda * this_min
// }

pub fn compute_data_costs<T, F>(
    parameters: &common::Parameters,
    data_cost_function: F,
) -> Vec<Vec<T>>
where
    F: Fn(u8, u8) -> T,
    // T: std::default::Default + Clone,
    T: std::default::Default + Copy,
{
    let mut data = Vec::with_capacity(parameters.total_pixels as usize);
    for _y in 0..parameters.height {
        for _x in 0..i64::from(parameters.maximum_disparity - 1) {
            data.push(vec![T::default(); parameters.maximum_disparity as usize]);
        }
        for x in i64::from(parameters.maximum_disparity - 1)..i64::from(parameters.width) {
            let mut current_pixel_data = Vec::with_capacity(parameters.maximum_disparity as usize);
            for d in 0..i64::from(parameters.maximum_disparity) {
                let data_cost = {
                    // if (x - d) < 0 {
                    //     data_cost_function(255u8, 0u8)
                    // } else {
                    //     data_cost_function(
                    //         parameters.left_image[x as usize],
                    //         parameters.right_image[(x - d) as usize],
                    //     )
                    // }
                    data_cost_function(
                        parameters.left_image[x as usize],
                        parameters.right_image[(x - d) as usize],
                    )
                };
                current_pixel_data.push(data_cost);
            }
            data.push(current_pixel_data);
        }
    }
    data
}

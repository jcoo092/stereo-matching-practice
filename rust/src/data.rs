use crate::common;
// use std::ops;

// pub fn squared_difference<T: ops::Mul<Output = T> + ops::Sub<Output = T> + Copy>(a: T, b: T) -> T {
//     (a - b) * (a - b)
//}

// pub fn absolute_difference<T: num_traits::sign::abs>(a: T, b: T) -> T {
//     a.abs(b)
// }

pub fn absolute_difference_u8_to_f32(a: u8, b: u8) -> f32 {
    (f32::from(a) - f32::from(b)).abs()
}

pub fn compute_data_costs<T, F>(
    parameters: &common::Parameters,
    data_cost_function: F,
) -> Vec<Vec<T>>
where
    F: Fn(u8, u8) -> T,
{
    let mut data = Vec::with_capacity(parameters.total_pixels as usize);
    for _y in 0..parameters.height {
        for x in i64::from(parameters.maximum_disparity - 1)..i64::from(parameters.width) {
            // let leftIdx = x + y * parameters.width;
            let mut current_pixel_data = Vec::with_capacity(parameters.maximum_disparity as usize);
            for d in 0..i64::from(parameters.maximum_disparity) {
                let data_cost = {
                    if (x - d) < 0 {
                        data_cost_function(255u8, 0u8)
                    } else {
                        data_cost_function(
                            parameters.left_image[x as usize],
                            parameters.right_image[(x - d) as usize],
                        )
                    }
                };
                current_pixel_data.push(data_cost);
            }
            data.push(current_pixel_data);
        }
    }
    data
}

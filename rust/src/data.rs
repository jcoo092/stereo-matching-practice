use crate::common;
use std::ops;

pub fn squared_difference<T: ops::Mul<Output = T> + ops::Sub<Output = T> + Copy>(a: T, b: T) -> T {
    (a - b) * (a - b)
}

// pub fn absolute_difference<T: num_traits::sign::abs>(a : T, b : T) -> T {
//     a.abs(b)
// }

pub fn absolute_difference_u8_to_f32(a: u8, b: u8) -> f32 {
    (a as f32 - b as f32).abs()
}

pub fn compute_data_costs<T, F>(
    parameters: &common::Parameters,
    data_cost_function: F,
) -> Vec<Vec<T>>
where
    F: Fn(u8, u8) -> T,
{
    let mut data = Vec::with_capacity(parameters.total_pixels as usize);
    for y in 0..parameters.height {
        for x in 0..parameters.width as i64 {
            // let leftIdx = x + y * parameters.width;
            let mut currentPixelData = Vec::with_capacity(parameters.maximum_disparity as usize);
            for d in 0..parameters.maximum_disparity as i64 {
                let data_cost = {
                    if (x - d) < 0 {
                        data_cost_function(255u8, 0u8)
                    } else {
                        data_cost_function(parameters.left_image[x as usize], parameters.right_image[(x - d) as usize])
                    }
                };
                currentPixelData.push(data_cost);
            }
            data.push(currentPixelData);
        }
    }
    data
}

use crate::common;

pub fn squared_difference<T>(a : T, b : T) -> T {
    (a - b) * (a - b)
}

pub fn compute_data_costs<T, F>(parameters : common::Parameters, data_cost_function : F) -> Vec<Vec<T>> where F: Fn(u8, u8) -> T {
    let mut data = Vec::with_capacity(parameters.total_pixels as usize);
    for y in 0..parameters.height {
        for x in 0..parameters.width {
            // let leftIdx = x + y * parameters.width;
            let mut currentPixelData = Vec::with_capacity(parameters.maximum_disparity as usize);
            for d in 0..parameters.maximum_disparity {
                let data_cost = {
                    if x - d < 0 {
                        data_cost_function 255 0
                    } else {
                        data_cost_function parameters.left_image[x] parameters.right_image[x - d]
                    }
                }
                currentPixelData.push(data_cost);
            }
            data.push(currentPixelData);
        }
    }
    data
}
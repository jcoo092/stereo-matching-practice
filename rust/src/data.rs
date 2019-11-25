use crate::common;

pub fn absolute_difference_u8_to_f32(a: u8, b: u8) -> f32 {
    (f32::from(a) - f32::from(b)).abs()
}

pub fn truncated_linear_f32_fh(lambda: f32, tau: f32, a: u8, b: u8) -> f32 {
    lambda * f32::min(tau, absolute_difference_u8_to_f32(a, b))
}

pub fn compute_data_costs<T, F>(
    parameters: &common::Parameters,
    data_cost_function: F,
) -> Vec<Vec<T>>
where
    F: Fn(u8, u8) -> T,
    T: Default + Copy,
{
    let mut data = Vec::with_capacity(parameters.total_pixels);
    for y in 0..parameters.height {
        for x in 0..(parameters.width) {
            let pixel_index = (x as usize) + (y as usize) * parameters.width;
            let mut current_pixel_data = Vec::with_capacity(parameters.maximum_disparity);
            for d in 0..(parameters.maximum_disparity) {
                let data_cost = {
                    if (x as i64 - d as i64) < 0 {
                        data_cost_function(255u8, 0u8)
                    } else {
                        data_cost_function(
                            parameters.left_image[pixel_index],
                            parameters.right_image[(pixel_index - d)],
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

#[cfg(test)]
mod tests {
    use super::truncated_linear_f32_fh;
    use crate::common;

    #[test]
    fn test_truncated_linear_two_ones() {
        let computation = truncated_linear_f32_fh(common::LAMBDA_FH, common::TAU_FH, 1, 1);
        let result = (computation - 0.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_one_zero() {
        let computation = truncated_linear_f32_fh(common::LAMBDA_FH, common::TAU_FH, 1, 0);
        let result = (computation - common::LAMBDA_FH).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_two_hundred_fifty_five_one() {
        let computation = truncated_linear_f32_fh(common::LAMBDA_FH, common::TAU_FH, 255, 1);
        let result = (computation - (common::LAMBDA_FH * common::TAU_FH)).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_two_hundred_fifty_five_zero() {
        let computation = truncated_linear_f32_fh(common::LAMBDA_FH, common::TAU_FH, 255, 0);
        let result = (computation - (common::LAMBDA_FH * common::TAU_FH)).abs();
        assert!(result <= std::f32::EPSILON);
    }
}

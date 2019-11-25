use crate::common;

pub fn truncated_linear_f32_fh(d: f32, a: usize, b: usize) -> f32 {
    f32::min((a as i32 - b as i32).abs() as f32, d)
}

pub fn compute_smoothness_costs<T, F>(
    parameters: &common::Parameters,
    smoothness_cost_function: F,
) -> Vec<Vec<T>>
where
    F: Fn(usize, usize) -> T,
{
    let mut smoothness_costs = Vec::with_capacity(parameters.maximum_disparity);
    for y in 0..parameters.maximum_disparity {
        let mut this_row = Vec::with_capacity(parameters.maximum_disparity);
        for x in 0..parameters.maximum_disparity {
            this_row.push(smoothness_cost_function(x, y));
        }
        smoothness_costs.push(this_row);
    }
    smoothness_costs
}

#[cfg(test)]
mod tests {
    use super::compute_smoothness_costs;
    use super::truncated_linear_f32_fh;
    use crate::common;

    #[test]
    fn test_truncated_linear_fh_identical() {
        let computation = truncated_linear_f32_fh(common::D_FH, 5, 5);
        let result = (computation - 0.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_one_apart_a() {
        let computation = truncated_linear_f32_fh(common::D_FH, 6, 5);
        let result = (computation - 1.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_one_apart_b() {
        let computation = truncated_linear_f32_fh(common::D_FH, 5, 6);
        let result = (computation - 1.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_two_apart_a() {
        let computation = truncated_linear_f32_fh(common::D_FH, 7, 5);
        let result = (computation - 1.7f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_two_apart_b() {
        let computation = truncated_linear_f32_fh(common::D_FH, 5, 7);
        let result = (computation - 1.7f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_ten_apart_a() {
        let computation = truncated_linear_f32_fh(common::D_FH, 11, 1);
        let result = (computation - 1.7f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_truncated_linear_fh_ten_apart_b() {
        let computation = truncated_linear_f32_fh(common::D_FH, 1, 11);
        let result = (computation - 1.7f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_compute_smoothness_costs_fh_five_by_five_d_two() {
        let params = common::Parameters {
            left_image: &Vec::new(),
            right_image: &Vec::new(),
            width: 5,
            height: 5,
            total_pixels: 25,
            window_edge_size: 3,
            maximum_disparity: 2,
            use_zero_mean: false,
        };

        let expected_result = vec![vec![0.0, 1.0], vec![1.0, 0.0]];

        let actual_result =
            compute_smoothness_costs(&params, |a, b| truncated_linear_f32_fh(common::D_FH, a, b));

        assert_eq!(actual_result, expected_result);
    }

    #[test]
    fn test_compute_smoothness_costs_fh_ten_by_ten_d_five() {
        let params = common::Parameters {
            left_image: &Vec::new(),
            right_image: &Vec::new(),
            width: 10,
            height: 10,
            total_pixels: 25,
            window_edge_size: 3,
            maximum_disparity: 5,
            use_zero_mean: false,
        };

        let expected_result = vec![
            vec![0.0, 1.0, 1.7, 1.7, 1.7],
            vec![1.0, 0.0, 1.0, 1.7, 1.7],
            vec![1.7, 1.0, 0.0, 1.0, 1.7],
            vec![1.7, 1.7, 1.0, 0.0, 1.0],
            vec![1.7, 1.7, 1.7, 1.0, 0.0],
        ];

        let actual_result =
            compute_smoothness_costs(&params, |a, b| truncated_linear_f32_fh(common::D_FH, a, b));

        assert_eq!(actual_result, expected_result);
    }
}

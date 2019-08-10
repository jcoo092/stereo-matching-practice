// Parameters common to every matching algorithm
pub struct Parameters<'a> {
    pub left_image: &'a Vec<u8>,
    pub right_image: &'a Vec<u8>,
    pub width: u32,
    pub height: u32,
    pub total_pixels: u32,
    pub window_edge_size: u32,
    pub maximum_disparity: u32,
    pub use_zero_mean: bool,
}

pub const LAMBDA_FH: f32 = 0.07;
pub const TAU_FH: f32 = 15.0;
pub const D_FH: f32 = 1.7;

// Blatantly stolen from StackOverflow...  unfortunately I have lost track of precisely where...
pub fn compute_mean_of_vec<T>(input_vec: &[T]) -> T
where
    for<'x> T: Copy
        + num::Zero
        + std::ops::Add<T, Output = T>
        + std::ops::Div<T, Output = T>
        + num::FromPrimitive
        + std::iter::Sum<&'x T>,
{
    let sum: T = input_vec.iter().sum();
    sum / num::FromPrimitive::from_usize(input_vec.len()).unwrap()
}

pub fn argmin_of_vec<T>(input_vec: &[T]) -> usize
where
    T: Copy + PartialOrd,
{
    let first_value = input_vec[0];
    let (result_index, _) = input_vec.iter().enumerate().fold(
        (0, first_value),
        |(index, smallest_so_far), (next_index, element)| {
            if *element < smallest_so_far {
                (next_index, *element)
            } else {
                (index, smallest_so_far)
            }
        },
    );
    result_index
}

#[cfg(test)]
mod tests {
    use super::argmin_of_vec;
    use super::compute_mean_of_vec;
    #[test]
    fn test_compute_mean_of_vec_f32_five_fives() {
        let test_vec = vec![5.0f32; 5];
        let mean = compute_mean_of_vec(&test_vec);
        let result = (mean - 5.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_compute_mean_of_vec_f32_one_to_five_ascending() {
        let test_vec: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5];
        let mean = compute_mean_of_vec(&test_vec);
        let result = (mean - 3.3f32).abs();
        assert!(result <= std::f32::EPSILON);
    }

    #[test]
    fn test_argmin_of_vec_f32_five_to_one_descending() {
        let test_vec: Vec<f32> = vec![5.5, 4.4, 3.3, 2.2, 1.1];
        let result = argmin_of_vec(&test_vec);
        assert_eq!(result, 4);
    }

    #[test]
    fn test_argmin_of_vec_f32_one_to_five_ascending() {
        let test_vec: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5];
        let result = argmin_of_vec(&test_vec);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_argmin_of_vec_f32_five_to_one_jumbled() {
        let test_vec: Vec<f32> = vec![3.3, 5.5, 1.1, 2.2, 4.4];
        let result = argmin_of_vec(&test_vec);
        assert_eq!(result, 2);
    }

    #[test]
    fn test_t_default_f32() {
        let actual_result: Vec<f32> = vec![Default::default(); 3];
        let expected_result = vec![0.0f32; 3];
        assert_eq!(actual_result, expected_result);
    }
}

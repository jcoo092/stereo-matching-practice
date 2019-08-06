// Parameters common to every matching algorithm
pub struct Parameters {
    pub left_image: Vec<u8>,
    pub right_image: Vec<u8>,
    pub width: u32,
    pub height: u32,
    pub total_pixels: u32,
    pub window_edge_size: u32,
    pub maximum_disparity: u32,
    pub use_zero_mean: bool,
}

// pub fn compute_mean_of_vec<
//     'a,
//     T: std::iter::Sum<&'a T> + std::convert::From<usize> + std::ops::Div<Output = T>,
// >(
//     input_vec: &'a Vec<T>,
// ) -> T {
//     input_vec.iter().sum::<T>() / T::from(input_vec.len())
// }

// pub fn compute_mean_of_vec<'a, T>(input_vec: &'a Vec<T>) -> T
// where
//     T: std::iter::Sum<&'a T> + std::convert::From<usize> + std::ops::Div<Output = T>,
// {
//     input_vec.iter().sum::<T>() / T::from(input_vec.len())
// }

// pub fn compute_mean_of_f32_vec(input_vec: &[f32]) -> f32 {
//     let total: f32 = input_vec.iter().sum();
//     total / (input_vec.len() as f32)
//     // input_vec.iter().sum() / (input_vec.len() as f32)
// }

// pub fn compute_mean_of_vec<T>(input_vec: &Vec<T>) -> T {
//     let total = input_vec.iter().sum();
//     total / T::from(input_vec.len() as u32)
// }

// Blatantly stolen from StackOverflow...  unfortunately I have lost track of precisely where...
// pub fn compute_mean_of_vec<'g, T>(input_vec: &'g Vec<T>) -> T
pub fn compute_mean_of_vec<T>(input_vec: &[T]) -> T
where
    for<'x> T: Copy
        + num::Zero
        + std::ops::Add<T, Output = T>
        + std::ops::Div<T, Output = T>
        + num::FromPrimitive
        + std::iter::Sum<&'x T>,
    // T: num::Num + num::FromPrimitive,
{
    let sum: T = input_vec.iter().sum();
    sum / num::FromPrimitive::from_usize(input_vec.len()).unwrap()
}

pub fn argmin_of_vec<T>(input_vec: &[T]) -> usize
where
    T: Copy + std::cmp::PartialOrd,
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
    // use super::compute_mean_of_f32_vec;
    use super::compute_mean_of_vec;
    #[test]
    // fn test_compute_mean_of_f32_vec() {
    //     let test_vec = vec![5.0f32, 5.0f32, 5.0f32, 5.0f32];
    //     let mean = compute_mean_of_f32_vec(&test_vec);
    //     let result = (mean - 5.0f32).abs();
    //     assert!(result <= std::f32::EPSILON);
    // }
    fn test_compute_mean_of_vec_f32_five_fives() {
        let test_vec = vec![5.0f32; 5];
        let mean = compute_mean_of_vec(&test_vec);
        let result = (mean - 5.0f32).abs();
        assert!(result <= std::f32::EPSILON);
    }
}

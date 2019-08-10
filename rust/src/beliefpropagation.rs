use crate::common;
use crate::data;
use crate::smoothness;
use std::convert::TryInto;

pub struct BPParameters<T> {
    pub number_of_iterations: u32,
    pub data_cost_function: fn(u8, u8) -> T,
    pub smoothness_cost_function: fn(u32, u32) -> T,
}

fn compute_neighbours(parameters: &common::Parameters, index: usize) -> Vec<usize> {
    let usize_width = parameters.width as usize;
    let x = index % usize_width;
    let y = index / usize_width;

    let mut neighbour_vec = Vec::new();
    if y > 0 {
        neighbour_vec.push(index - usize_width);
    }
    if x > 0 {
        neighbour_vec.push(index - 1);
    }
    if x < (usize_width - 1) {
        neighbour_vec.push(index + 1);
    }
    if y < (parameters.height as usize - 1) {
        neighbour_vec.push(index + usize_width);
    }

    neighbour_vec
}

fn normalise_cost_vec<T>(cost_vec: &mut Vec<T>)
where
    for<'x> T: Copy
        + num::traits::identities::Zero
        + num::traits::cast::FromPrimitive
        + std::iter::Sum<&'x T>
        + num::traits::NumOps
        + num::traits::NumAssignOps,
{
    let mean = common::compute_mean_of_vec(cost_vec);
    for c in cost_vec.iter_mut() {
        *c -= mean;
    }
}

// fn normalise_cost_vec<T>(cost_vec: &mut Vec<T>)
// where
//     for<'x> T:
//         Copy + num::traits::NumAssignOps + std::cmp::PartialOrd + num::traits::identities::Zero,
// {
//     let min = common::min_partial_ord(cost_vec);
//     if min != T::zero() {
//         for c in cost_vec.iter_mut() {
//             *c /= min;
//         }
//     }
// }

fn normalise_all_messages<'b, T>(messages: &'b mut Vec<Vec<Vec<T>>>)
where
    for<'x> T: Copy
        + num::traits::cast::FromPrimitive
        + num::traits::identities::Zero
        + std::iter::Sum<&'x T>
        + num::traits::NumOps
        + num::traits::NumAssignOps
        + std::cmp::PartialOrd,
{
    for m in messages.iter_mut() {
        for w in m.iter_mut() {
            normalise_cost_vec(w);
        }
    }
}

// #[allow(dead_code)]
// #[cfg(test)]
fn find_index_in_neighbour(this_index: usize, neighbourhood: &[usize]) -> usize {
    neighbourhood
        .iter()
        .position(|x| *x == this_index)
        .expect("Didn't find the current pixel in its neighbour's neighbourhood")
}

fn update_messages<
    T: std::default::Default
        + std::clone::Clone
        + Copy
        + std::cmp::PartialOrd
        + num::traits::NumOps
        + num::traits::NumAssignOps,
>(
    parameters: &common::Parameters,
    data_costs: &[Vec<T>],
    smoothness_costs: &[Vec<T>],
    neighbourhoods: &[Vec<usize>],
    messages1: &[Vec<Vec<T>>],
    messages2: &mut Vec<Vec<Vec<T>>>,
) {
    let max_d = parameters.maximum_disparity as usize;
    for (i, m1) in messages1.iter().enumerate() {
        // let neighbour_message_sums_start = std::time::Instant::now();
        let mut neighbour_messages_sums = vec![T::default(); max_d];
        for (fp, nms) in neighbour_messages_sums.iter_mut().enumerate() {
            for m1ni in m1.iter() {
                *nms += data_costs[i][fp] + m1ni[fp];
            }
        }
        // let neighbour_message_sums_end = std::time::Instant::now();
        // println!(
        //     "data costs computation time took {:?}.",
        //     neighbour_message_sums_end.duration_since(neighbour_message_sums_start)
        // );

        for (j, neighbour_index) in neighbourhoods[i].iter().enumerate() {
            // let index_in_neighbour_start = std::time::Instant::now();
            // let index_in_neighbour = neighbourhoods[*neighbour_index]
            //     .iter()
            //     .position(|x| *x == i)
            //     .expect("Didn't find the current pixel in its neighbour's neighbourhood");
            let index_in_neighbour = find_index_in_neighbour(i, &neighbourhoods[*neighbour_index]);
            // let index_in_neighbour_end = std::time::Instant::now();
            // println!(
            //     "index in neighbour computation time took {:?}.",
            //     index_in_neighbour_end.duration_since(index_in_neighbour_start)
            // );

            let m2_neighbour_costs = &mut messages2[*neighbour_index][index_in_neighbour];
            for (fq, m2nc) in m2_neighbour_costs.iter_mut().enumerate() {
                let min_cost = (0..max_d)
                    .map(|fp| smoothness_costs[fp][fq] + neighbour_messages_sums[fp] - m1[j][fp])
                    .min_by(|a, b| a.partial_cmp(b).expect("Found a NaN in this vector..."))
                    .unwrap(); // If this is ever empty, something has REALLY gone wrong
                *m2nc = min_cost;
            }
        }
    }
}

fn compute_final_disparity<T>(
    max_d: usize,
    data_costs: &[Vec<T>],
    pixel_index: usize,
    pixel_messages: &[Vec<T>],
) -> u8
where
    T: std::default::Default + std::ops::AddAssign + Copy + std::cmp::PartialOrd,
{
    let mut beliefs: Vec<T> = Vec::with_capacity(max_d);
    for i in 0..max_d {
        let mut message_cost = data_costs[pixel_index][i];
        for neighbour_costs in pixel_messages {
            message_cost += neighbour_costs[i];
        }
        beliefs.push(message_cost);
    }

    common::argmin_of_vec(&beliefs).try_into().unwrap()
}

pub fn belief_propagation<'c, T>(
    parameters: &common::Parameters,
    bpparameters: &BPParameters<T>,
) -> Vec<u8>
where
    for<'x> T: 'c
        + std::default::Default
        + std::clone::Clone
        + Copy
        + std::cmp::PartialOrd
        + num::traits::cast::FromPrimitive
        + num::traits::identities::Zero
        + std::iter::Sum<&'x T>
        + num::traits::NumOps
        + num::traits::NumAssignOps
        + std::fmt::Debug,
{
    // let data_costs_start = std::time::Instant::now();
    let data_costs = data::compute_data_costs(&parameters, bpparameters.data_cost_function);
    // let data_costs_end = std::time::Instant::now();
    // println!(
    //     "data costs computation time took {:?}.",
    //     data_costs_end.duration_since(data_costs_start)
    // );

    // let smoothness_costs_start = std::time::Instant::now();
    let smoothness_costs =
        smoothness::compute_smoothness_costs(&parameters, bpparameters.smoothness_cost_function);
    // let smoothness_costs_end = std::time::Instant::now();
    // println!(
    //     "smoothness costs computation time took {:?}.",
    //     smoothness_costs_end.duration_since(smoothness_costs_start)
    // );

    // let neighbourhoods_start = std::time::Instant::now();
    let neighbourhoods = {
        let mut neighbour_vecs = vec![Vec::<usize>::new(); parameters.total_pixels as usize];
        for (i, v) in neighbour_vecs.iter_mut().enumerate() {
            *v = compute_neighbours(&parameters, i)
        }
        neighbour_vecs
    };
    // let neighbourhoods_end = std::time::Instant::now();
    // println!(
    //     "neighbourhoods computation time took {:?}.",
    //     neighbourhoods_end.duration_since(neighbourhoods_start)
    // );

    // let create_messages_start = std::time::Instant::now();
    let mut messages1 = (0..parameters.total_pixels as usize)
        .map(|i| {
            vec![vec![T::default(); parameters.maximum_disparity as usize]; neighbourhoods[i].len()]
        })
        .collect::<Vec<_>>();

    let mut messages2 = messages1.clone();
    // let create_messages_end = std::time::Instant::now();
    // println!(
    //     "messages creation time took {:?}.",
    //     create_messages_end.duration_since(create_messages_start)
    // );

    for _i in 0..bpparameters.number_of_iterations {
        // let update_messages_start = std::time::Instant::now();
        update_messages(
            parameters,
            &data_costs,
            &smoothness_costs,
            &neighbourhoods,
            &messages1,
            &mut messages2,
        );
        // let update_messages_end = std::time::Instant::now();
        // println!(
        //     "update messages computation time took {:?}.",
        //     update_messages_end.duration_since(update_messages_start)
        // );

        // for m in messages2.iter().skip(60000).take(2) {
        //     println!("messages2 pre-normalisation: {:?}", m);
        // }

        // let normalise_messages_start = std::time::Instant::now();
        normalise_all_messages(&mut messages2);
        // let normalise_messages_end = std::time::Instant::now();
        // println!(
        //     "normalise messages computation time took {:?}.",
        //     normalise_messages_end.duration_since(normalise_messages_start)
        // );

        // for m in messages2.iter().skip(60000).take(2) {
        //     println!("messages2 post-normalisation: {:?}", m);
        // }

        // let swap_messages_start = std::time::Instant::now();
        std::mem::swap(&mut messages1, &mut messages2);
        // messages2 = std::mem::replace(&mut messages1, messages2);
        // let swap_messages_end = std::time::Instant::now();
        // println!(
        //     "normalise messages computation time took {:?}.",
        //     swap_messages_end.duration_since(swap_messages_start)
        // );
    }

    // let final_disparities_start = std::time::Instant::now();
    messages1
        .iter()
        .enumerate()
        .map(|(i, p)| {
            compute_final_disparity(parameters.maximum_disparity as usize, &data_costs, i, p)
        })
        .collect()
    // let final_disparities_end = std::time::Instant::now();
    // println!(
    //     "final disparities computation time took {:?}.",
    //     final_disparities_end.duration_since(final_disparities_start)
    // );
}

#[cfg(test)]
mod tests {
    use super::compute_neighbours;
    use super::find_index_in_neighbour;
    use super::normalise_cost_vec;
    use crate::common;
    // use std::convert::TryInto;

    // #[test]
    // fn test_compute_final_disparity_a() {
    //     let sample_costs = [[
    //         [
    //             -1.548_749_9,
    //             -0.5687499,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //             0.15125012,
    //         ],
    //         [
    //             -1.5499998,
    //             -0.54999983,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.150_000_21,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //             0.15000021,
    //         ],
    //     ]];
    //     let actual_result = sample_costs
    //     .iter()
    //     .enumerate()
    //     .map(|(i, p)| {
    //         compute_final_disparity(parameters.maximum_disparity as usize, &data_costs, i, p)
    //     })
    //     .collect()
    // }

    #[test]
    fn test_compute_neighbours_four_rows_three_columns_a() {
        let params = common::Parameters {
            left_image: &Vec::new(),
            right_image: &Vec::new(),
            width: 3,
            height: 4,
            total_pixels: 4 * 3,
            window_edge_size: 3,
            maximum_disparity: 2,
            use_zero_mean: false,
        };

        let actual_neighbourhoods = (0..params.total_pixels as usize)
            .map(|i| compute_neighbours(&params, i))
            .collect::<Vec<_>>();

        let expected_neighbourhoods = vec![
            vec![1, 3],
            vec![0, 2, 4],
            vec![1, 5],
            vec![0, 4, 6],
            vec![1, 3, 5, 7],
            vec![2, 4, 8],
            vec![3, 7, 9],
            vec![4, 6, 8, 10],
            vec![5, 7, 11],
            vec![6, 10],
            vec![7, 9, 11],
            vec![8, 10],
        ];

        assert_eq!(actual_neighbourhoods, expected_neighbourhoods);
    }

    #[test]
    fn test_compute_neighbours_four_rows_three_columns_b() {
        let params = common::Parameters {
            left_image: &Vec::new(),
            right_image: &Vec::new(),
            width: 3,
            height: 4,
            total_pixels: 4 * 3,
            window_edge_size: 3,
            maximum_disparity: 2,
            use_zero_mean: false,
        };

        let actual_neighbourhoods = {
            let mut neighbour_vecs = vec![Vec::<usize>::new(); params.total_pixels as usize];
            for (i, v) in neighbour_vecs.iter_mut().enumerate() {
                *v = compute_neighbours(&params, i)
            }
            neighbour_vecs
        };

        let expected_neighbourhoods = vec![
            vec![1, 3],
            vec![0, 2, 4],
            vec![1, 5],
            vec![0, 4, 6],
            vec![1, 3, 5, 7],
            vec![2, 4, 8],
            vec![3, 7, 9],
            vec![4, 6, 8, 10],
            vec![5, 7, 11],
            vec![6, 10],
            vec![7, 9, 11],
            vec![8, 10],
        ];

        assert_eq!(actual_neighbourhoods, expected_neighbourhoods);
    }

    // #[test]
    // fn test_normalise_cost_vec_f32_one_to_five_ascending() {
    //     let mut test_vec: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5];
    //     normalise_cost_vec(&mut test_vec);
    //     let expected_result: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5]
    //         .iter()
    //         .map(|x| *x - 3.3)
    //         .collect();
    //     // let result = (mean - 3.3f32).abs();
    //     assert_eq!(test_vec, expected_result);
    // }

    #[test]
    fn test_normalise_cost_vec_f32_one_to_five_ascending() {
        let mut test_vec: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5];
        normalise_cost_vec(&mut test_vec);
        let expected_result: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5]
            .iter()
            .map(|x| *x - 3.3)
            .collect();
        // let result = (mean - 3.3f32).abs();
        assert_eq!(test_vec, expected_result);
    }

    #[test]
    fn test_find_index_in_neighbour_four_rows_three_columns() {
        let params = common::Parameters {
            left_image: &Vec::new(),
            right_image: &Vec::new(),
            width: 3,
            height: 4,
            total_pixels: 4 * 3,
            window_edge_size: 3,
            maximum_disparity: 2,
            use_zero_mean: false,
        };

        let actual_neighbourhoods = (0..params.total_pixels as usize)
            .map(|i| compute_neighbours(&params, i))
            .collect::<Vec<_>>();

        let actual_indices: Vec<Vec<usize>> = actual_neighbourhoods
            .iter()
            .enumerate()
            .map(|(i, neighbourhood)| {
                neighbourhood
                    .iter()
                    .map(|j| find_index_in_neighbour(i, &actual_neighbourhoods[*j]))
                    .collect()
            })
            .collect();

        let expected_indices = vec![
            vec![0, 0],
            vec![0, 0, 0],
            vec![1, 0],
            vec![1, 1, 0],
            vec![2, 1, 1, 0],
            vec![1, 2, 0],
            vec![2, 1, 0],
            vec![3, 1, 1, 0],
            vec![2, 2, 0],
            vec![2, 1],
            vec![3, 1, 1],
            vec![2, 2],
        ];

        assert_eq!(actual_indices, expected_indices);
    }
}

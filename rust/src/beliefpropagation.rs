use crate::common;
use crate::data;
use crate::smoothness;
use std::convert::TryInto;

pub struct BPParameters<T> {
    pub number_of_iterations: usize,
    pub data_cost_function: fn(u8, u8) -> T,
    pub smoothness_cost_function: fn(usize, usize) -> T,
}

fn compute_neighbours(parameters: &common::Parameters, index: usize) -> Vec<usize> {
    let usize_width = parameters.width;
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
    if y < (parameters.height - 1) {
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

fn normalise_all_messages<'b, T>(messages: &'b mut Vec<Vec<Vec<T>>>)
where
    for<'x> T: Copy
        + num::traits::cast::FromPrimitive
        + num::traits::identities::Zero
        + std::iter::Sum<&'x T>
        + num::traits::NumOps
        + num::traits::NumAssignOps
        + PartialOrd,
{
    for m in messages.iter_mut() {
        for w in m.iter_mut() {
            normalise_cost_vec(w);
        }
    }
}

fn find_index_in_neighbour(this_index: usize, neighbourhood: &[usize]) -> usize {
    neighbourhood
        .iter()
        .position(|x| *x == this_index)
        .expect("Didn't find the current pixel in its neighbour's neighbourhood")
}

fn update_messages<
    T: Default + Clone + Copy + PartialOrd + num::traits::NumOps + num::traits::NumAssignOps,
>(
    parameters: &common::Parameters,
    data_costs: &[Vec<T>],
    smoothness_costs: &[Vec<T>],
    neighbourhoods: &[Vec<usize>],
    messages1: &[Vec<Vec<T>>],
    messages2: &mut Vec<Vec<Vec<T>>>,
) {
    let max_d = parameters.maximum_disparity;
    for (i, m1) in messages1.iter().enumerate() {
        let mut neighbour_messages_sums = vec![T::default(); max_d];
        for (fp, nms) in neighbour_messages_sums.iter_mut().enumerate() {
            for m1ni in m1.iter() {
                *nms += data_costs[i][fp] + m1ni[fp];
            }
        }

        for (j, neighbour_index) in neighbourhoods[i].iter().enumerate() {
            let index_in_neighbour = find_index_in_neighbour(i, &neighbourhoods[*neighbour_index]);

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
    T: Default + std::ops::AddAssign + Copy + PartialOrd,
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
        + Default
        + Copy
        + PartialOrd
        + num::traits::cast::FromPrimitive
        + num::traits::identities::Zero
        + std::iter::Sum<&'x T>
        + num::traits::NumOps
        + num::traits::NumAssignOps
        + std::fmt::Debug,
{
    let data_costs = data::compute_data_costs(&parameters, bpparameters.data_cost_function);

    let smoothness_costs =
        smoothness::compute_smoothness_costs(&parameters, bpparameters.smoothness_cost_function);

    let neighbourhoods = {
        let mut neighbour_vecs = vec![Vec::<usize>::new(); parameters.total_pixels];
        for (i, v) in neighbour_vecs.iter_mut().enumerate() {
            *v = compute_neighbours(&parameters, i)
        }
        neighbour_vecs
    };

    let mut messages1 = (0..parameters.total_pixels)
        .map(|i| {
            vec![vec![T::default(); parameters.maximum_disparity]; neighbourhoods[i].len()]
        })
        .collect::<Vec<_>>();

    let mut messages2: Vec<Vec<Vec<T>>> = messages1.clone();

    for _i in 0..bpparameters.number_of_iterations {
        update_messages(
            parameters,
            &data_costs,
            &smoothness_costs,
            &neighbourhoods,
            &messages1,
            &mut messages2,
        );

        normalise_all_messages(&mut messages2);

        std::mem::swap(&mut messages1, &mut messages2);
    }

    messages1
        .iter()
        .enumerate()
        .map(|(i, p)| {
            compute_final_disparity(parameters.maximum_disparity, &data_costs, i, p)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::compute_neighbours;
    use super::find_index_in_neighbour;
    use super::normalise_cost_vec;
    use crate::common;

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

        let actual_neighbourhoods = (0..params.total_pixels)
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
            let mut neighbour_vecs = vec![Vec::<usize>::new(); params.total_pixels];
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

    #[test]
    fn test_normalise_cost_vec_f32_one_to_five_ascending() {
        let mut test_vec: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5];
        normalise_cost_vec(&mut test_vec);
        let expected_result: Vec<f32> = vec![1.1, 2.2, 3.3, 4.4, 5.5]
            .iter()
            .map(|x| *x - 3.3)
            .collect();
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

        let actual_neighbourhoods = (0..params.total_pixels)
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

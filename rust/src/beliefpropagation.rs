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

// fn normalise_cost_f32_vec(cost_vec: &mut Vec<f32>) {
//     let average = common::compute_mean_of_f32_vec(cost_vec);
//     for c in cost_vec {
//         *c -= average;
//     }
// }

fn normalise_cost_vec<T>(cost_vec: &mut Vec<T>)
where
    for<'x> T: std::ops::SubAssign
        + Copy
        + num::traits::identities::Zero
        + std::ops::Div<Output = T>
        + num::traits::cast::FromPrimitive
        + std::iter::Sum<&'x T>,
{
    let mean = common::compute_mean_of_vec(cost_vec);
    for c in cost_vec.iter_mut() {
        *c -= mean;
    }
}

// fn normalise_all_messages_f32(message_vec: &mut Vec<Vec<Vec<f32>>>) {
//     for v in message_vec.iter_mut() {
//         for w in v.iter_mut() {
//             normalise_cost_f32_vec(w);
//         }
//     }
// }

fn normalise_all_messages<'b, T>(message_vec: &'b mut Vec<Vec<Vec<T>>>)
where
    for<'x> T: Copy
        + std::ops::SubAssign
        + num::traits::cast::FromPrimitive
        + std::ops::Div<Output = T>
        + num::traits::identities::Zero
        + std::iter::Sum<&'x T>,
{
    for v in message_vec.iter_mut() {
        for w in v.iter_mut() {
            normalise_cost_vec(w);
        }
    }
}

fn update_messages<
    T: std::default::Default
        + std::ops::Add<Output = T>
        + std::clone::Clone
        + Copy
        + std::ops::AddAssign
        + std::ops::Sub<Output = T>
        + std::cmp::PartialOrd,
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
        let mut neighbour_messages_sums = vec![T::default(); max_d];
        for (fp, nms) in neighbour_messages_sums.iter_mut().enumerate() {
            // for neighbourhoods_index in 0..neighbourhoods[i].len() {
            for m1ni in m1.iter() {
                *nms += data_costs[i][fp] + m1ni[fp];
            }
        }
        for (j, neighbour_index) in neighbourhoods[i].iter().enumerate() {
            // let neighbour_costs = &m1[*neighbour_index];
            // let index_in_neighbour = [*neighbour_index];
            // let neighbour_costs = &m1[index_in_neighbour];
            let index_in_neighbour = neighbourhoods[*neighbour_index]
                .iter()
                .position(|x| *x == i)
                .expect("Didn't find the current pixel in its neighbour's neighbourhood");
            let m2_neighbour_costs = &mut messages2[*neighbour_index][index_in_neighbour];
            // for fq in 0..max_d {
            //     // let mutable min_cost =
            //     let min_cost = (0..max_d)
            //         .map(|fp| {
            //             smoothness_costs[fp][fq] + neighbour_messages_sums[fp] - neighbour_costs[fp]
            //         })
            //         .min_by(|a, b| a.partial_cmp(b).expect("Found a NaN in this vector..."));
            //     m2_neighbour_costs[fq] = min_cost.unwrap(); // If this is ever empty, something has REALLY gone wrong
            // }
            for (fq, m2nc) in m2_neighbour_costs.iter_mut().enumerate() {
                let min_cost = (0..max_d)
                    .map(|fp| smoothness_costs[fp][fq] + neighbour_messages_sums[fp] - m1[j][fp])
                    .min_by(|a, b| a.partial_cmp(b).expect("Found a NaN in this vector..."))
                    .unwrap();
                // m2_neighbour_costs[fq] = min_cost.unwrap(); // If this is ever empty, something has REALLY gone wrong
                *m2nc = min_cost; // If this is ever empty, something has REALLY gone wrong
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
        // let data_cost = &data_costs[pixel_index][i];
        // let mut message_cost = T::default();
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
        + std::ops::Add<Output = T>
        + Copy
        + std::ops::Sub<Output = T>
        + std::ops::AddAssign
        + std::cmp::PartialOrd
        + std::ops::SubAssign
        + num::traits::cast::FromPrimitive
        + num::traits::identities::Zero
        + std::ops::Div<Output = T>
        + std::iter::Sum<&'x T>,
{
    let data_costs = data::compute_data_costs(&parameters, bpparameters.data_cost_function);
    let smoothness_costs =
        smoothness::compute_smoothness_costs(&parameters, bpparameters.smoothness_cost_function);

    let neighbourhoods = {
        let mut neighbour_vecs = vec![Vec::<usize>::new(); parameters.total_pixels as usize];
        for (i, v) in neighbour_vecs.iter_mut().enumerate() {
            *v = compute_neighbours(&parameters, i)
        }
        neighbour_vecs
    };

    let mut messages1 = (0..parameters.total_pixels as usize)
        .map(|i| {
            let this_pixels_neighbours = &neighbourhoods[i];
            this_pixels_neighbours
                .iter()
                .map(|_| vec![T::default(); parameters.maximum_disparity as usize])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut messages2 = messages1.clone();

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
            compute_final_disparity(parameters.maximum_disparity as usize, &data_costs, i, p)
        })
        .collect()
}

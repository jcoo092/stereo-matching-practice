use crate::common;
use crate::data;
use crate::smoothness;

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

fn normalise_cost_f32_vec(cost_vec: &mut Vec<f32>) {
    let average = common::compute_mean_of_f32_vec(cost_vec);
    for c in cost_vec {
        *c -= average;
    }
}

fn normalise_cost_vec<'a, T: 'a>(cost_vec: &'a mut Vec<T>) where T: std::ops::SubAssign + Copy + num::traits::identities::Zero + std::ops::Div<Output = T> + num::traits::cast::FromPrimitive + std::iter::Sum<&'a T> {
    let average = common::compute_mean_of_vec(cost_vec);
    for c in cost_vec {
        *c -= average;
    }
}

fn normalise_all_messages_f32(message_vec: &mut Vec<Vec<Vec<f32>>>) {
    for v in message_vec.iter_mut() {
        for w in v.iter_mut() {
            normalise_cost_f32_vec(w);
        }
    }
}

fn normalise_all_messages<'a, T: 'a>(message_vec: &'a mut Vec<Vec<Vec<T>>>) where T: Copy + std::ops::SubAssign + num::traits::cast::FromPrimitive + std::ops::Div<Output = T> + num::traits::identities::Zero + std::iter::Sum<&'a T> {
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
    for (i, v) in messages1.iter().enumerate() {
        let mut neighbour_messages_sums = vec![T::default(); max_d];
        for (fp, nms) in neighbour_messages_sums.iter_mut().enumerate() {
            for neighbour_index in &neighbourhoods[i] {
                // neighbour_messages_sums[fp] =
                // neighbour_messages_sums[fp] + data_costs[i][fp] + v[*neighbour_index][fp];
                *nms += data_costs[i][fp] + v[*neighbour_index][fp];
            }
        }
        for (index_in_neighbour, neighbour_index) in neighbourhoods[i].iter().enumerate() {
            let neighbour_costs = &v[*neighbour_index];
            let m2_neighbour_costs = &mut messages2[*neighbour_index][index_in_neighbour];
            for fq in 0..max_d {
                // let mutable min_cost =
                let min_cost = (0..max_d)
                    .map(|fp| {
                        smoothness_costs[fp][fq] + neighbour_messages_sums[fp] - neighbour_costs[fp]
                    })
                    .min_by(|a, b| a.partial_cmp(b).expect("Found a NaN in this vector..."));
                m2_neighbour_costs[fq] = min_cost.unwrap(); // If this is ever empty, something has REALLY gone wrong
            }
        }
    }
}

pub fn belief_propagation<'a,
    T: 'a +  std::default::Default
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
        + std::iter::Sum<&'a T>
>(
    parameters: &common::Parameters,
    bpparameters: &BPParameters<T>,
) -> Vec<u8> {
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

    // let mut messages1 = vec![
    //     vec![T::default(); parameters.maximum_disparity as usize];
    //     parameters.total_pixels as usize
    // ];

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

    Vec::new()
}

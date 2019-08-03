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

fn update_messages<
    T: std::default::Default + std::ops::Add<Output = T> + std::clone::Clone + Copy,
>(
    parameters: &common::Parameters,
    bpparameters: &BPParameters<T>,
    data_costs: &[Vec<T>],
    smoothness_costs: &[Vec<T>],
    neighbourhoods: &[Vec<usize>],
    messages1: &[Vec<T>],
    messages2: &mut Vec<Vec<T>>,
) {
    let maxD = parameters.maximum_disparity as usize;
    for (i, v) in messages1.iter().enumerate() {
        let mut neighbour_messages_sums = vec![T::default(); maxD];
        for fp in 0..maxD {
            for neighbour_index in &neighbourhoods[i] {
                neighbour_messages_sums[fp] = neighbour_messages_sums[fp]
                    + data_costs[i][fp]
                    + messages1[*neighbour_index][fp];
            }
        }
    }
}

pub fn belief_propagation<
    T: std::default::Default + std::clone::Clone + std::ops::Add<Output = T> + Copy,
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

    let mut messages1 = vec![
        vec![T::default(); parameters.maximum_disparity as usize];
        parameters.total_pixels as usize
    ];

    let mut messages2 = messages1.clone();

    for _i in 0..bpparameters.number_of_iterations {
        update_messages(
            parameters,
            bpparameters,
            &data_costs,
            &smoothness_costs,
            &neighbourhoods,
            &messages1,
            &mut messages2,
        );

        std::mem::swap(&mut messages1, &mut messages2);
    }

    Vec::new()
}

use crate::common;
use crate::data;

pub struct BPParameters<T> {
    pub number_of_iterations: u32,
    pub data_cost_function: fn(u8, u8) -> T,
    smoothness_function: fn(u32, u32) -> T
}

pub fn belief_propagation<T>(parameters: &common::Parameters, bpparameters: &BPParameters<T>) -> Vec<u8> {
    let data_costs = data::compute_data_costs(&parameters, bpparameters.data_cost_function);
    let smoothness_costs = smoothness::compute_smoothness_costs(&parameters, bpparameters.smoothness_function);

    Vec::new()
}

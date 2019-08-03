use crate::common;

pub const D_FH : f32 = 1.7;

pub fn truncated_linear(d: f32, a: u32, b: u32) -> f32 {
    f32::min((a as i32 - b as i32).abs() as f32, d)
}

pub fn compute_smoothness_costs<T, F>(parameters: &common::Parameters, smoothness_cost_function: F) -> Vec<Vec<T>> where F: Fn(u32, u32) -> T, {
    let mut smoothness_costs = Vec::with_capacity(parameters.maximum_disparity as usize);
    for x in 0..parameters.maximum_disparity {
        let mut this_row = Vec::with_capacity(parameters.maximum_disparity as usize);
        for y in 0..parameters.maximum_disparity {
            this_row.push(smoothness_cost_function(x, y));
        }
        smoothness_costs.push(this_row);
    }
    smoothness_costs
}
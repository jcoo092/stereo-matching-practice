use crate::common;
use crate::data;

pub struct BPParameters {
    pub number_of_iterations: u32,
}

pub fn belief_propagation(parameters: &common::Parameters, bpparameters: &BPParameters) -> Vec<u8> {
    // unimplemented!();
    parameters.left_image.to_owned()
}

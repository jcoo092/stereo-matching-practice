// Parameters common to every matching algorithm
pub struct Parameters {
    pub left_image: Vec<u8>,
    pub right_image: Vec<u8>,
    pub width: u32,
    pub height: u32,
    pub total_pixels: u32,
    pub window_edge_size: usize,
    pub maximum_disparity: usize,
    pub use_zero_mean: bool,
}

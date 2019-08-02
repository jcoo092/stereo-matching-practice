use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug)]
enum Algorithms {
    SAD,
    SSD,
    DynamicProgramming,
    BeliefPropagation,
}

impl std::str::FromStr for Algorithms {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        //let algorithm =
        match s.to_lowercase().as_str() {
            "sad" => Ok(Algorithms::SAD),
            "ssd" => Ok(Algorithms::SSD),
            "dynamicprogramming" | "dynamic_programming" | "dynamic-programming" | "dp" => {
                Ok(Algorithms::DynamicProgramming)
            }
            "beliefpropagation" | "belief_propagation" | "belief-propagation" | "bp" => {
                Ok(Algorithms::BeliefPropagation)
            }
            _ => Err(
                "Parsing of the algorithm failed.  Please specify one of the possible options."
                    .to_string(),
            ),
        }

        //Ok(algorithm)
    }
}

impl std::fmt::Display for Algorithms {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let algo_string = match self {
            Algorithms::SAD => "sad",
            Algorithms::SSD => "ssd",
            Algorithms::DynamicProgramming => "dynamic-programming",
            Algorithms::BeliefPropagation => "belief-propagation",
        };
        write!(f, "{}", algo_string)
    }
}

#[derive(Debug, StructOpt)]
#[structopt(
    name = "Practice stereo matcher",
    about = "A command-line program to execute various stereo matching algorithms on stereo pairs of images."
)]
struct CLIParameters {
    #[structopt(
        short = "l",
        parse(from_os_str),
        help = "Path to the left image to be matched on."
    )]
    left_image_path: PathBuf,
    #[structopt(
        short = "r",
        parse(from_os_str),
        help = "Path to the right image to be matched on."
    )]
    right_image_path: PathBuf,
    #[structopt(
        short = "o",
        parse(from_os_str),
        help = "Path to the output directory for the disparity map to be output to.\nNote that you DO NOT provide a filename, just the directory.  The filename is generated by the program itself."
    )]
    output_directory: PathBuf,
    #[structopt(
        short = "w",
        help = "The size of a side of the square window you would like to use for the matching.  Currently only relevant to SAD and SSD."
    )]
    window_size: Option<usize>,
    #[structopt(
        short = "d",
        help = "The maximum disparity size that the program will search out to (defaults to 32)."
    )]
    maximum_disparity: Option<usize>,
    #[structopt(
        short = "a",
        help = "The algorithm you would like to use for the stereo matching."
    )]
    algorithm: Algorithms,
    #[structopt(
        short = "n",
        help = "The number of iterations of the algorithm to carry out.  Currently only relevant to belief propagation."
    )]
    number_of_iterations: Option<usize>,
    #[structopt(
        short = "z",
        help = "Set this flag if you would like the program to use the zero-mean version of the unary cost function.  (Currently does nothing)"
    )]
    use_zero_mean: bool,
}

fn determine_output_file_path(params: &CLIParameters) -> PathBuf {
    let algorithm_string = std::ffi::OsString::from(format!("_{}", params.algorithm));
    let left_image_name_without_extension = params
        .left_image_path
        .file_stem()
        .expect("Left image name apparently has no file name (???)");
    let left_image_extension = params
        .left_image_path
        .extension()
        .expect("Left image name apparently has no extension (???)");
    let window_size = {
        let base_string: String = if params.window_size.is_some() {
            format!("_{}", params.window_size.unwrap())
        } else {
            String::new()
        };
        std::ffi::OsString::from(base_string)
    };

    let number_of_iterations = {
        let base_string = if params.number_of_iterations.is_some() {
            format!("_{}", params.number_of_iterations.unwrap())
        } else {
            String::new()
        };
        std::ffi::OsString::from(base_string)
    };

    [
        left_image_name_without_extension,
        &algorithm_string,
        &window_size,
        &number_of_iterations,
        left_image_extension,
    ]
    .iter()
    .collect()
}

fn main() {
    let cli_parameters = CLIParameters::from_args();
    println!("{:?}", cli_parameters);
    println!("{:#?}", cli_parameters.window_size);
    //println!("Hello, world!");
}

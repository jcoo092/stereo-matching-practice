// Learn more about F# at http://fsharp.org

open System
open Argu

type MatchingAlgorithms =
    | SAD
    | SSD
    | DynamicProgramming
    | BeliefPropagation

type CLIArguments =
    | [<Mandatory;AltCommandLine("-l")>]LeftImage of left_image_path:string
    | [<Mandatory;AltCommandLine("-r")>]RightImage of right_image_path:string
    | [<Mandatory;AltCommandLine("-o")>]OutputDirectory of output_directory:string
    | [<EqualsAssignment;AltCommandLine("-w")>]Window of window_size:int
    | [<EqualsAssignment;AltCommandLine("-d")>]MaximumDisparity of max_disparity:int
    | [<Mandatory;AltCommandLine("-a");>]Algorithm of matching_algorithm:MatchingAlgorithms
    | [<CliPrefix(CliPrefix.Dash)>] Z
with
    interface IArgParserTemplate with
        member s.Usage =
            match s with
            | LeftImage _ -> "Path (including filename) to the left image to be matched"
            | RightImage _ -> "Path (including filename) to the right image to be matched"
            | OutputDirectory _ -> "The directory for the file to be output to.
Note that you do NOT specify a filename - that will be automatically constructed based on the filename of the left image and the options selected"
            | Window _ -> "Size of the sides of the square window to be used, assuming that the selected algorithm uses windows.
The total size of the window will thus be n x n, where n is the input size here"
            | MaximumDisparity _ -> "The maximum disparity to search across, assuming that the selected algorithm uses a maximum"
            | Algorithm _ -> "The specific choice of stereo matching algorithm to be used"
            | Z -> "Use the zero-mean version of the algorithm (only applies currently to SAD and SSD)"

let processWindowSizeArgument (results : ParseResults<CLIArguments>) =
    let result = results.TryGetResult Window
    match result with
    | Some(w) -> w
    | None -> 3

let processMaximumDisparityArgument (results : ParseResults<CLIArguments>) =
    let result = results.TryGetResult Window
    match result with
    | Some(d) -> d
    | None -> 32

let doNothing = 5

[<EntryPoint>]
let main argv =

    let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some ConsoleColor.Red)
    let parser = ArgumentParser.Create<CLIArguments>(programName = "stereo matching runner", errorHandler = errorHandler)

    let results = parser.ParseCommandLine argv

    printfn "Got parse results %A" <| results.GetAllResults()

    let returnedArray =
        match results.GetResult Algorithm with
        | SAD -> doNothing
        | SSD -> doNothing
        | DynamicProgramming -> doNothing
        | BeliefPropagation -> doNothing





    //printfn "Hello World from F#!"
    0 // return an integer exit code

// Learn more about F# at http://fsharp.org

open System
open Argu
open System.IO
open Common
open SAD
open DynamicProgramming
open SixLabors.ImageSharp
open SixLabors.ImageSharp.PixelFormats
open SixLabors.ImageSharp.Processing
open SixLabors.ImageSharp.Advanced

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

let getAlgorithmString =
    function
    | SAD -> "SAD"
    | SSD -> "SSD"
    | DynamicProgramming -> "Dynamic_Programming"
    | BeliefPropagation -> "Belief_Propagation"

let determineOutputFilename (results : ParseResults<CLIArguments>) =
    let algorithmString = results.GetResult Algorithm |> getAlgorithmString
    let leftImageName = results.GetResult LeftImage
    let leftimagewithoutextension = Path.GetFileNameWithoutExtension leftImageName
    let leftImageExtension = Path.GetExtension leftImageName
    let windowSize = if results.Contains Window then "_" + (results.GetResult Window |> string) else String.Empty
    //leftimagewithoutextension + "_" + algorithmString + leftImageExtension
    sprintf "%s_%s%s.%s" leftimagewithoutextension algorithmString windowSize leftImageExtension

let openImageAndConvertToGrayscaleArray (imagePath : string) =
    use img = Image.Load(imagePath)
    img.Mutate(fun x -> x.Grayscale() |> ignore)
    //let pspan = img.GetPixelSpan()
    img.GetPixelSpan().ToArray() |> Array.Parallel.map (fun p -> p.R)
    //img.GetPixelSpan()

let getImageSize (imagePath : string) =
    use img = Image.Load(imagePath)
    (img.Width, img.Height)

let makeGray8 intensity = Gray8(intensity)

[<EntryPoint>]
let main argv =

    let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some ConsoleColor.Red)
    let parser = ArgumentParser.Create<CLIArguments>(programName = "stereo matching runner", errorHandler = errorHandler)

    let results = parser.ParseCommandLine argv

    //printfn "Got parse results %A" <| results.GetAllResults()

    let imgWidth, imgHeight = results.GetResult LeftImage |> getImageSize

    let outputImageArray =
        let matchingParameters = {
            leftImage = results.GetResult LeftImage |> openImageAndConvertToGrayscaleArray
            rightImage = results.GetResult RightImage |> openImageAndConvertToGrayscaleArray
            width = imgWidth
            height = imgHeight
            windowEdgeSize = results.TryGetResult Window |> Option.defaultValue 3
            maximumDisparity = results.TryGetResult MaximumDisparity |> Option.defaultValue 32
            zeroMean = results.Contains Z
        }
        match results.GetResult Algorithm with
        | SAD -> raise (NotImplementedException "This stereo matching algorithm has not yet been implemented")
        | SSD -> raise (NotImplementedException "This stereo matching algorithm has not yet been implemented")
        | DynamicProgramming -> 
            let updatedMatchingParameters = {
                leftImage = matchingParameters.leftImage |> Array.Parallel.map uint32
                rightImage = matchingParameters.rightImage |> Array.Parallel.map uint32
                width = matchingParameters.width
                height = matchingParameters.height
                windowEdgeSize = matchingParameters.windowEdgeSize
                maximumDisparity = matchingParameters.maximumDisparity
                zeroMean = matchingParameters.zeroMean
            }
            dynamicProgramming updatedMatchingParameters |> Array.Parallel.map byte
        | BeliefPropagation -> raise (NotImplementedException "This stereo matching algorithm has not yet been implemented")

    let outputImage = Image.LoadPixelData(Array.Parallel.map makeGray8 outputImageArray, imgWidth, imgHeight)

    let outputFilename = (results.GetResult OutputDirectory) + (string Path.DirectorySeparatorChar) +
                            (determineOutputFilename results)

    //use outFile = new System.IO.FileStream((results.GetResult OutputDirectory) + @"\" + outputFilename, FileMode.OpenOrCreate)

    outputImage.Save(outputFilename)

    printfn "Saved stereo-matching-result image to %s" outputFilename
    //printfn "Hello World from F#!"
    0 // return an integer exit code

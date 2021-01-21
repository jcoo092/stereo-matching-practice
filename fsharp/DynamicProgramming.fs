// NB:  This was never finished (so far).  It does NOT work, and should not be referred to as an example of effective 1D Dynamic Programming Stereo Matching
// If that is a topic you are interested in, I strongly recommend searching online for papers and other materials by (now Emeritus Professor) Georgy Gimel'farb
// who first published on this particular topic back in the 1970's, and is quite possibly one of the world's top experts on the topic.
module DynamicProgramming

open Common
open System

let diffWholeLine (leftLine: 'a []) (rightLine: 'a []) x d =
    // printfn "line length is %d, x is %d and d is %d.  x - d = %d" (Array.length rightLine) x d (x - d)

    Data.carefulAbsoluteDifference leftLine.[x] rightLine.[x - d]
    |> uint32

let matchAlongScanline maxDisparity differenceFunction (leftLine: 'a [], rightLine: 'a []) =
    let diffF = differenceFunction leftLine rightLine
    let lineLength = leftLine.Length

    let forwardPass =
        Array2D.zeroCreate lineLength (maxDisparity + 1)

    let backwardPass =
        Array2D.zeroCreate lineLength (maxDisparity + 1)

    let endOfLine = lineLength - 1
    let mutable emin = UInt32.MaxValue
    let mutable dmin = 0u

    forwardPass.[0, 0] <- diffF 0 0

    for x = 1 to endOfLine do
        dmin <- 0u

        // for d = 0 to ((min (x - 1) maxDisparity) - 1) do
        for d = 0 to ((min (x - 1) maxDisparity)) do
            // emin <- UInt32.MaxValue

            // for d' = (max 0 (d - 1)) to ((min (x - 2) (maxDisparity))) do
            //     // for d' = (max 0 (d - 1)) to ((min (x - 2) (maxDisparity)) - 1) do
            //     let error = forwardPass.[(x - 1), d']
            //     // let (p, q) = minAndArgmin forwardPass.[(x - 1), *]

            //     if error < emin then
            //         emin <- error
            //         dmin <- uint32 <| d'

            let (em, dm) =
                minAndArgmin forwardPass.[(x - 1), (max 0 (d - 1))..(min (x - 2) (maxDisparity))]

            // forwardPass.[x, d] <- emin + (diffF x d)
            // backwardPass.[x, d] <- dmin

            forwardPass.[x, d] <- em + (diffF x d)
            backwardPass.[x, d] <- uint32 <| dm

    emin <- UInt32.MaxValue
    dmin <- 0u

    for d in 0 .. ((min endOfLine maxDisparity)) do
        // for d in 0 .. ((min endOfLine maxDisparity) - 1) do
        if forwardPass.[endOfLine, d] < emin then
            emin <- forwardPass.[endOfLine, d]
            // dmin <- d |> uint32
            dmin <- uint32 <| d


    let outputArray = Array.zeroCreate lineLength
    outputArray.[endOfLine] <- dmin

    for x = (endOfLine - 1) downto 0 do
        // for x = (endOfLine - 1) downto 0 do
        // let opi = outputArray.[x + 1] |> int
        let opi = int <| outputArray.[x + 1]
        outputArray.[x] <- backwardPass.[x + 1, opi]

    outputArray

let dynamicProgramming parameters =
    let leftSlices = buildArraySlices parameters Left
    let rightSlices = buildArraySlices parameters Right
    let slicesZip = Array.zip leftSlices rightSlices

    let computed =
        Array.map (fun l -> matchAlongScanline parameters.maximumDisparity diffWholeLine l) slicesZip
    // Array.Parallel.map (fun l -> matchAlongScanline parameters.maximumDisparity diffWholeLine l) slicesZip

    Array.concat computed

// NB:  This was never finished (so far).  It does NOT work, and should not be referred to as an example of effective 1D Dynamic Programming Stereo Matching
// If that is a topic you are interested in, I strongly recommend searching online for papers and other materials by (now Emeritus Professor) Georgy Gimel'farb
// who first published on this particular topic back in the 1970's, and is quite possibly one of the world's top experts on the topic.
module DynamicProgramming

open Common
open System

let matchAlongScanline maxDisparity differenceFunction (leftLine : 'a [], rightLine : 'a []) =
    let diffF = differenceFunction leftLine rightLine
    let lineLength = leftLine.Length
    let forwardPass = Array2D.zeroCreate (maxDisparity + 1) lineLength
    let backwardPass = Array2D.zeroCreate (maxDisparity + 1) lineLength

    let endOfLine = lineLength - 1
    let mutable emin = UInt32.MaxValue
    let mutable dmin = 0u

    forwardPass.[0,0] <- diffF 0 0
    for x = 1 to endOfLine do
        for d = 0 to (min (x - 1) maxDisparity) do
            emin <- UInt32.MaxValue
            for d' = (max 0 (d - 1)) to (min (x - 2) (maxDisparity)) do
                let error = diffF (x - 1) d'
                if error < emin then
                    emin <- error
                    dmin <- d' |> uint32

            forwardPass.[x,d] <- emin + diffF x d
            backwardPass.[x, d] <- dmin

    emin <- UInt32.MaxValue
    dmin <- 0u

    for d in 0..(min endOfLine maxDisparity) do
        if forwardPass.[endOfLine, d] < emin then
            emin <- forwardPass.[endOfLine, d]
            dmin <- d |> uint32

    let outputArray = Array.zeroCreate lineLength
    outputArray.[endOfLine] <- dmin

    for x = (endOfLine - 1) downto 0 do
        let opi = outputArray.[x + 1] |> int
        outputArray.[x] <- backwardPass.[x + 1, opi]

    outputArray

let dynamicProgramming parameters =
    let leftSlices = buildArraySlices parameters Left
    let rightSlices = buildArraySlices parameters Right
    let _slicesZip = Array.zip leftSlices rightSlices
    [||]
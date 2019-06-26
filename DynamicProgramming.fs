module DynamicProgramming

open Common
open System

let matchAlongScanline maxDisparity differenceFunction (leftLine : 'a [], rightLine : 'a []) =
    let diffF = differenceFunction leftLine rightLine
    let lineLength = leftLine.Length
    let forwardPass = Array2D.zeroCreate (maxDisparity + 1) lineLength
    let backwardPass = Array2D.zeroCreate (maxDisparity + 1) lineLength

    let endOfLine = lineLength - 1
    //let mutable emin = Byte.MaxValue
    //let mutable dmin = 0uy
    let mutable emin = UInt32.MaxValue
    let mutable dmin = 0u

    forwardPass.[0,0] <- diffF 0 0
    //for x in 1..endOfLine do
    for x = 1 to endOfLine do
        //for d in 0..(min (x - 1) maxDisparity) do
        for d = 0 to (min (x - 1) maxDisparity) do
            //emin <- Byte.MaxValue
            emin <- UInt32.MaxValue
            //let lb = max 0 (d - 1)
            //let ub = min (x - 2) (maxDisparity)
            //for d' in lb..ub do
            for d' = (max 0 (d - 1)) to (min (x - 2) (maxDisparity)) do
                let error = diffF (x - 1) d'
                if error < emin then
                    emin <- error
                    //dmin <- byte d'
                    dmin <- d' |> uint32

            forwardPass.[x,d] <- emin + diffF x d
            backwardPass.[x, d] <- dmin

    //emin <- Byte.MaxValue
    //dmin <- 0uy
    emin <- UInt32.MaxValue
    dmin <- 0u

    for d in 0..(min endOfLine maxDisparity) do
        if forwardPass.[endOfLine, d] < emin then
            emin <- forwardPass.[endOfLine, d]
            //dmin <- byte d
            dmin <- d |> uint32

    let outputArray = Array.zeroCreate lineLength
    outputArray.[endOfLine] <- dmin

    for x = (endOfLine - 1) downto 0 do
        let opi = outputArray.[x + 1] |> int
        outputArray.[x] <- backwardPass.[x + 1, opi]

    outputArray

let dynamicProgramming parameters =
    //let leftSlices = buildSlices parameters Left
    //let leftSlices = buildArraySegments parameters Left
    let leftSlices = buildArraySlices parameters Left
    //let rightSlices = buildSlices parameters Right
    //let rightSlices = buildArraySegments parameters Right
    let rightSlices = buildArraySlices parameters Right
    let slicesZip = Array.zip leftSlices rightSlices
    //let matchedLines = Array.Parallel.map (matchAlongScanline parameters.maximumDisparity Data.arraysSquaredDifference) slicesZip
    //let matchedLines = Array.Parallel.map (matchAlongScanline parameters.maximumDisparity slicesSquaredDifference) slicesZip
    //Array.concat matchedLines
    [||]
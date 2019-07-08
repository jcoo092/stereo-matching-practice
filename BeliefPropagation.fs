// Currently just trying to focus on a basic version of min-sum 2D loopy belief propagation
// Using min-sum on the same basis that Felzenswalb & Huttenlocher substituted max-product with it

module BeliefPropagation

open Common
open System

//[<Struct; IsByRefLike>]
type BPParameters = {
    dataFunction: Byte -> Byte -> float32
    smoothnessFunction: int -> int -> float32
    iterations: int
}

[<Struct>]
type Direction =
    | North = 0
    | East = 1
    | South = 2
    | West = 3
    | Data = 4

type Pixel<'a> = {
    msg : 'a[,]
    mutable bestDisparity: int
}

// type MRF2D<'a> = {
//     grid : Pixel<'a>[]
//     width: int
//     height: int
// }

// let computeNeighbours parameters i =
//     let x = i % parameters.width
//     let y = i / parameters.width

//     [|
//         if y > 0 then yield i - parameters.width;
//         if x > 0 then yield i - 1;
//         if x < (parameters.width - 1) then yield i + 1;
//         if y < (parameters.height - 1) then yield i + parameters.width;
//     |]

// let inline initMessages parameters =
//     Array.Parallel.init (parameters.width * parameters.height) (
//         fun i ->
//                 let neighbours = computeNeighbours parameters i
//                 Array.map (fun neighbour -> (neighbour, Array.zeroCreate (parameters.maximumDisparity + 1))) neighbours
//     )

let initMessages parameters (dataCosts : 'a[][]) =
    Array.Parallel.init (parameters.width * parameters.height) (
        fun i ->
            let message = Array2D.create 5 parameters.maximumDisparity 0.0f
            //Array.blit dataCosts.[i] 0 message.[int Direction.Data, *] 0 parameters.maximumDisparity
            for j = 0 to (parameters.maximumDisparity - 1) do
                message.[int Direction.Data, j] <- dataCosts.[i].[j]
            {msg = message; bestDisparity = 0}
    )

// let inline normalizeCostArray arr =
//     let sum = Array.min arr
//     Array.iteri (fun i value -> arr.[i] <- value / sum) arr

// // This is intended to match eq. 2 in F & H 2006
// let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) (m1 : (int * float32 []) [] []) (m2 : (int * float32 []) [] []) =
// // this function computes updates to the messages using data in m1, and stores it back to m2
// // p and q are used below in accordance with Felzenswalb & Huttenlocher's notation
//     Array.Parallel.iteri (fun i p -> // each pixel in the image
//         let fpMax = min maxD ((Array.length dataCosts.[i]) - 1)
//         // let fpMax = ((Array.length dataCosts.[i]) - 1)
//         let neighbourMessageSums = Array.zeroCreate (fpMax + 1)
//         for fp = 0 to fpMax do
//             for (_, (neighbourCosts : float32 [])) in p do
//                 neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]
//                 // Strictly speaking, data costs shouldn't be here, but since it is all additions, and data costs vary only by fp
//                 // It's easier just to include them here
//         for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
//             let indexInNeighbour : int = Array.findIndex (fst >> ((=) i)) m1.[neighbourIdx]
//             let (_, m2neighbourcosts) = m2.[neighbourIdx].[indexInNeighbour]
//             for fq = 0 to maxD do // each disparity label of q
//                 let mutable mincost = Single.MaxValue
//                 for fp = 0 to fpMax do
//                     let smoothnessCost = smoothnessCosts.[fp, fq]
//                     let previousMessageCost = neighbourMessageSums.[fp] - neighbourCosts.[fp]
//                     let totalCost = smoothnessCost + previousMessageCost
//                     if totalCost < mincost then
//                         mincost <- totalCost


//                 m2neighbourcosts.[fq] <- mincost
//     ) m1
//     Array.Parallel.iter (fun p ->
//                         Array.iter (fun (_, neighbourCosts) -> normalizeCostArray neighbourCosts) p
//     ) m2

// let computeFinalDisparities parameters (dataCosts : float32 [][]) (messages : (int * float32 []) [] []) =
//     Array.Parallel.mapi (fun i p ->
//         let maxFq = min parameters.maximumDisparity ((Array.length dataCosts.[i]) - 1)

//         // Compute belief vector
//         let beliefs = Array.zeroCreate (maxFq + 1)
//         for j = 0 to maxFq do
//             let dataCost = dataCosts.[i].[j]
//             let mutable messageCost = 0.0f
//             for (_, (neighbourArray : float32 [])) in p do
//                 messageCost <- messageCost + neighbourArray.[j]
//             beliefs.[j] <- dataCost + messageCost

//         // Select disparity value with minimum cost
//         argminFloat32Array beliefs |> byte
//     ) messages

// let computeEnergy (dataCosts : float32 [][]) (smoothnessCosts : float32[,]) (messages : (int * float32 []) [] []) (finalDisparities : byte[]) =
//     let dC = Array.fold (fun acc i ->
//                             let finDepI = finalDisparities.[i] |> int
//                             acc + dataCosts.[i].[finDepI]
//                         ) 0.0f [|0..(Array.length finalDisparities) - 1|]
//     let sC = Array.Parallel.mapi (fun i p ->
//                             let fp = finalDisparities.[i] |> int
//                             let mutable totalCost = 0.0f
//                             for (neighbourIdx, _) in p do
//                                 let fq = finalDisparities.[neighbourIdx] |> int
//                                 totalCost <- totalCost + smoothnessCosts.[fp, fq]
//                             totalCost
//                         ) messages |> Array.sum
//     dC + sC

let MAP parameters (smoothnessCosts : float32[,]) (pixels : Pixel<float32>[]) =
    let width = parameters.width
    let height = parameters.height
    Array.Parallel.iteri
        (fun i (x : Pixel<float32>) ->
            let mutable best = Single.MaxValue
            for j = 0 to (parameters.maximumDisparity - 1) do
                let cost = Array.sum x.msg.[*,j]
                if cost < best then
                    best <- cost
                    x.bestDisparity <- j
        ) pixels

    let mutable energy = 0.0f

    for y = 0 to (height - 1) do
        for x = 0 to (width - 1) do
            //let currentPixel = pixels.[x + y * width]
            let currentBest = pixels.[x + y * width].bestDisparity
            energy <- energy + pixels.[x + y * width].msg.[int Direction.Data, currentBest]

            if x - 1 >= 0 then
                energy <- energy + smoothnessCosts.[currentBest, pixels.[(x - 1) + y * width].bestDisparity]
            if x + 1 < width then
                energy <- energy + smoothnessCosts.[currentBest, pixels.[(x + 1) + y * width].bestDisparity]
            if y - 1 >= 0 then
                energy <- energy + smoothnessCosts.[currentBest, pixels.[x + (y - 1) * width].bestDisparity]
            if y + 1 < height then
                energy <- energy + smoothnessCosts.[currentBest, pixels.[x + (y + 1) * width].bestDisparity]

    energy


let sendMsg (parameters : Common.Parameters) (smoothnessCosts : float32[,]) (messages : Pixel<float32>[]) x y (direction : Direction) =
    //let newMsg = Array.zeroCreate parameters.maximumDisparity
    let newMsg = Array2D.create 1 parameters.maximumDisparity 0.0f
    let width = parameters.width
    for i = 0 to (parameters.maximumDisparity - 1) do
        let mutable minVal = Single.MaxValue
        for j = 0 to (parameters.maximumDisparity - 1) do
            let mutable p = smoothnessCosts.[i, j] + messages.[x + y * width].msg.[int Direction.Data, j]
            //let mutable p = 0.0f
            //p <- p + smoothnessCosts.[i, j]
            //p <- p + messages.[x + y * width].msg.[int Direction.Data, j]

            if direction <> Direction.West then
                p <- p + messages.[x + y * width].msg.[int Direction.West, j]
            if direction <> Direction.East then
                p <- p + messages.[x + y * width].msg.[int Direction.East, j]
            if direction <> Direction.North then
                p <- p + messages.[x + y * width].msg.[int Direction.North, j]
            if direction <> Direction.South then
                p <- p + messages.[x + y * width].msg.[int Direction.South, j]

            // if (x + y * width) < 2 && i < 2  && j < 3 then
            //     printfn "pixNum = %d, i = %d, j = %d, minVal = %f, p = %f" (x + y * width) i j minVal p
            //     printfn "smoothness: %f, data: %f" smoothnessCosts.[i,j] messages.[x + y * width].msg.[int Direction.Data, j]
                //printfn "Data cost = %f" messages.[x + y * width].msg.[int Direction.Data, j]




            minVal <- min minVal p

        //newMsg.[i] <- minVal
        newMsg.[0,i] <- minVal

    // normalize newMsg
    // let minVal = Array.min newMsg
    // Array.Parallel.iteri (fun i x -> newMsg.[i] <- x / minVal) newMsg
    //let minVal = Array.min newMsg.[0,*]
    let minVal =
        let mini = Array.min newMsg.[0,*]
        if mini < 1.0f then
            1.0f
        else
            mini
    Array2D.iteri (fun i j x -> newMsg.[i,j] <- x / minVal) newMsg

    match direction with
    | Direction.West ->
        //Array.blit newMsg 0 messages.[y*width + (x - 1)].msg.[int Direction.East, *] 0 parameters.maximumDisparity
        Array2D.blit newMsg 0 0 messages.[y*width + (x - 1)].msg (int Direction.East) 0 1 parameters.maximumDisparity
    | Direction.East ->
        //Array.blit newMsg 0 messages.[y*width + (x + 1)].msg.[int Direction.West, *] 0 parameters.maximumDisparity
        Array2D.blit newMsg 0 0 messages.[y*width + (x + 1)].msg (int Direction.West) 0 1 parameters.maximumDisparity
    | Direction.North ->
        //Array.blit newMsg 0 messages.[(y - 1)*width + x].msg.[int Direction.South, *] 0 parameters.maximumDisparity
        Array2D.blit newMsg 0 0 messages.[(y - 1)*width + x].msg (int Direction.South) 0 1 parameters.maximumDisparity
    | Direction.South ->
        //Array.blit newMsg 0 messages.[(y + 1)*width + x].msg.[int Direction.North, *] 0 parameters.maximumDisparity
        Array2D.blit newMsg 0 0 messages.[(y + 1)*width + x].msg (int Direction.North) 0 1 parameters.maximumDisparity
    | _ -> ()

let bp (parameters : Common.Parameters) smoothnessCosts messages direction =
    match direction with
    | Direction.East ->
        for y = 0 to (parameters.height - 1) do
            for x = 0 to (parameters.width - 2) do
                sendMsg parameters smoothnessCosts messages x y direction
    | Direction.West ->
        for y = 0 to (parameters.height - 1) do
            for x = (parameters.width - 1) downto 1 do
                sendMsg parameters smoothnessCosts messages x y direction
    | Direction.South ->
        for x = 0 to (parameters.width - 1) do
            for y = 0 to (parameters.height - 2) do
                sendMsg parameters smoothnessCosts messages x y direction
    | Direction.North ->
        for x = 0 to (parameters.width - 1) do
            for y = (parameters.height - 1) downto 1 do
                sendMsg parameters smoothnessCosts messages x y direction
    | _ -> ()


let beliefpropagation parameters bpparameters =
    let dataCosts = Data.computeDataCosts parameters bpparameters.dataFunction
    let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters.smoothnessFunction
    let messages = initMessages parameters dataCosts // All messages will be 0 initially
    for i = 1 to bpparameters.iterations do
        // updateMessages parameters.maximumDisparity dataCosts smoothnessCosts messages1 messages2
        // let temp = messages1
        // messages1 <- messages2
        // messages2 <- temp
        bp parameters smoothnessCosts messages Direction.East
        bp parameters smoothnessCosts messages Direction.West
        bp parameters smoothnessCosts messages Direction.North
        bp parameters smoothnessCosts messages Direction.South

        let energy = MAP parameters smoothnessCosts messages

        printfn "iteration %d, energy = %f" i energy

    //let findeps = computeFinalDisparities parameters dataCosts messages
    //let finenergy = computeEnergy dataCosts smoothnessCosts messages1 findeps
    // printfn "Final energy is: %f" (finenergy / float32 parameters.totalPixels)
    //findeps

    Array.map (fun m -> m.bestDisparity |> byte) messages
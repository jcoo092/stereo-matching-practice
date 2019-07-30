// Currently just trying to focus on a basic version of min-sum 2D loopy belief propagation
// Using min-sum on the same basis that Felzenswalb & Huttenlocher substituted max-product with it

module BeliefPropagation

open Common
open System

[<Struct; Runtime.CompilerServices.IsByRefLike>]
type BPParameters = {
    dataFunction: Byte -> Byte -> float32
    smoothnessFunction: int -> int -> float32
    iterations: int
}

let hFH (dataCosts : float32 [][]) (messages : (int * float32 []) [] []) i exclusionIndex =
    let dcs = dataCosts.[i]
    let thisPixel = messages.[i]
    let neighbourCosts = Array.zeroCreate (Array.length dcs)
    for (_, neighbour) in thisPixel do
        Array.iteri (fun i x -> neighbourCosts.[i] <- neighbourCosts.[i] + x) neighbour
    Array.iteri (fun i x -> neighbourCosts.[i] <- x + dcs.[i] - (snd thisPixel.[exclusionIndex]).[i]) neighbourCosts
    neighbourCosts


// Function to compute the message update value using the Potts model
// Derived from the description in section 3.1 of F&H 2006
// let PottsFH (dataCosts : float32 [][]) (messages : (int * float32 []) [] []) i d minhfp fq =
//     let thisH = hFH dataCosts messages i

//     min (hFH dataCosts messages i fq) (minhfp + d)

let computeNeighbours parameters i =
    let x = i % parameters.width
    let y = i / parameters.width

    [|
        if y > 0 then yield i - parameters.width;
        if x > 0 then yield i - 1;
        if x < (parameters.width - 1) then yield i + 1;
        if y < (parameters.height - 1) then yield i + parameters.width;
    |]

let inline initMessages parameters (dataCosts : float32 [][]) =
    //Array.Parallel.init (parameters.width * parameters.height) (
    Array.init (parameters.width * parameters.height) (
        fun i ->
                let neighbours = computeNeighbours parameters i
                let arraysize = min ((parameters.maximumDisparity + 1)) ((Array.length dataCosts.[i]))
                Array.map (fun neighbour -> (neighbour, Array.zeroCreate arraysize)) neighbours
    )

let inline normalizeCostArray arr =
    let normalizer = Array.average arr
    Array.iteri (fun i value -> arr.[i] <- value - normalizer) arr

// This is pretty much directly lifted from F&H, both the paper and their sample code
let dt m =
    for fq = 1 to ((Array.length m) - 1) do
        m.[fq] <- min m.[fq] (m.[fq - 1] + Smoothness.C_FH)

    for fq = ((Array.length m) - 2) downto 0 do
        m.[fq] <- min m.[fq] (m.[fq + 1] + Smoothness.C_FH)


// This is intended to match eq. 2 in F & H 2006
//let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) (m1 : (int * float32 []) [] []) (m2 : (int * float32 []) [] []) =
let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) oddOrEven (m : (int * float32 []) [] []) =
// p and q are used below in accordance with Felzenswalb & Huttenlocher's notation
    // let startIndex = if oddOrEven then 0 else 1
    // //Array.Parallel.iteri (fun i p -> // each pixel in the image
    // for i in startIndex..2..((Array.length m) - 1) do
    //     let p = m.[i]
    // //Array.iteri (fun i p -> // each pixel in the image
    //     let fpMax = min maxD ((Array.length dataCosts.[i]) - 1)
    //     // let fpMax = ((Array.length dataCosts.[i]) - 1)
    //     let neighbourMessageSums = Array.zeroCreate (fpMax + 1)
    //     for fp = 0 to fpMax do
    //         for (_, (neighbourCosts : float32 [])) in p do
    //             neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]
    //             // Strictly speaking, data costs shouldn't be here, but since it is all additions, and data costs vary only by fp
    //             // It's easier just to include them here
    //     for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
    //         let indexInNeighbour : int = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
    //         let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
    //         for fq = 0 to min maxD ((Array.length mneighbourcosts) - 1) do // each disparity label of q
    //             let mutable mincost = Single.MaxValue
    //             for fp = 0 to fpMax do
    //                 let smoothnessCost = smoothnessCosts.[fp, fq]
    //                 let previousMessageCost = neighbourMessageSums.[fp] - neighbourCosts.[fp]
    //                 let totalCost = smoothnessCost + previousMessageCost
    //                 if totalCost < mincost then
    //                     mincost <- totalCost


    //             mneighbourcosts.[fq] <- mincost

    // Using F&H's 'min convolution' style, where the message from p to q is computed as
    // min_{f_p} (V(f_p - f_q) + h(f_p)), where h(f_p) is
    // D_p(f_p) + the sum of the messages for f_p received from the other neighbours at t-1
    let startIndex = if oddOrEven then 0 else 1
    //Array.Parallel.iteri (fun i p -> // each pixel in the image
    for i in startIndex..2..((Array.length m) - 1) do
        let p = m.[i]
    //Array.iteri (fun i p -> // each pixel in the image
        let fpMax = min maxD ((Array.length dataCosts.[i]) - 1)
        // let fpMax = ((Array.length dataCosts.[i]) - 1)
        let neighbourMessageSums = Array.zeroCreate (fpMax + 1)
        for fp = 0 to fpMax do
            for (_, (neighbourCosts : float32 [])) in p do
                neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]
                // Strictly speaking, data costs shouldn't be here, but since it is all additions, and data costs vary only by fp
                // It's easier just to include them here
        for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
            // printfn "i = %d, neighbourIdx = %d"  i neighbourIdx
            let hfp = Array.map2 (( - )) neighbourMessageSums neighbourCosts
            let minhfp = (Array.min hfp) + Smoothness.D_FH
            dt hfp
            let indexInNeighbour = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
            let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
            // printfn "mneighbourcosts length = %d, hfp length = %d, fpMax = %d" (Array.length mneighbourcosts) (Array.length hfp) fpMax
            for fq = 0 to (min fpMax (Array.length mneighbourcosts)) - 1 do // each disparity label of q
                // let mutable mincost = Single.MaxValue
                // for fp = 0 to fpMax do
                //     let smoothnessCost = smoothnessCosts.[fp, fq]
                //     let previousMessageCost = neighbourMessageSums.[fp] - neighbourCosts.[fp]
                //     let totalCost = smoothnessCost + previousMessageCost
                //     if totalCost < mincost then
                //         mincost <- totalCost
                mneighbourcosts.[fq] <- min hfp.[fq] minhfp



    //     // F&H Potts (sect 3.1)
    // let startIndex = if oddOrEven then 0 else 1
    // //Array.Parallel.iteri (fun i p -> // each pixel in the image
    // for i in startIndex..2..((Array.length m) - 1) do
    //     let p = m.[i]
    // //Array.iteri (fun i p -> // each pixel in the image
    //     let fpMax = min maxD ((Array.length dataCosts.[i]) - 1)
    //     // let fpMax = ((Array.length dataCosts.[i]) - 1)
    //     let neighbourMessageSums = Array.zeroCreate (fpMax + 1)
    //     for fp = 0 to fpMax do
    //         for (_, (neighbourCosts : float32 [])) in p do
    //             neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]
    //             // Strictly speaking, data costs shouldn't be here, but since it is all additions, and data costs vary only by fp
    //             // It's easier just to include them here
    //     for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
    //         // printfn "i = %d, neighbourIdx = %d"  i neighbourIdx
    //         let hfp = Array.map2 (( - )) neighbourMessageSums neighbourCosts
    //         let minhfp = (Array.min hfp) + Smoothness.D_FH
    //         let indexInNeighbour = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
    //         let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
    //         // printfn "mneighbourcosts length = %d, hfp length = %d, fpMax = %d" (Array.length mneighbourcosts) (Array.length hfp) fpMax
    //         for fq = 0 to (min fpMax (Array.length mneighbourcosts)) - 1 do // each disparity label of q
    //             // let mutable mincost = Single.MaxValue
    //             // for fp = 0 to fpMax do
    //             //     let smoothnessCost = smoothnessCosts.[fp, fq]
    //             //     let previousMessageCost = neighbourMessageSums.[fp] - neighbourCosts.[fp]
    //             //     let totalCost = smoothnessCost + previousMessageCost
    //             //     if totalCost < mincost then
    //             //         mincost <- totalCost
    //             mneighbourcosts.[fq] <- min hfp.[fq] minhfp


    //) m
    //Array.Parallel.iter (fun p ->
    Array.iter (fun p ->
                        Array.iter (fun (_, neighbourCosts) -> normalizeCostArray neighbourCosts) p
    ) m

let computeFinalDisparities parameters (dataCosts : float32 [][]) (messages : (int * float32 []) [] []) =
    // Array.Parallel.mapi (fun i p ->
    Array.mapi (fun i p ->
        let maxFq = min parameters.maximumDisparity ((Array.length dataCosts.[i]) - 1)

        // Compute belief vector
        let beliefs = Array.zeroCreate (maxFq + 1)
        for j = 0 to maxFq do
            let dataCost = dataCosts.[i].[j]
            let mutable messageCost = 0.0f
            for (_, (neighbourArray : float32 [])) in p do
                messageCost <- messageCost + neighbourArray.[j]
            beliefs.[j] <- dataCost + messageCost

        // Select disparity value with minimum cost
        argminFloat32Array beliefs |> byte
    ) messages

// let computeEnergy (dataCosts : float32 [][]) (smoothnessCosts : float32[,]) (messages : (int * float32 []) [] []) (finalDisparities : byte[]) =
//     let dC = Array.fold (fun acc i ->
//                             let finDepI = finalDisparities.[i] |> int
//                             acc + dataCosts.[i].[finDepI]
//                         ) 0.0f [|0..(Array.length finalDisparities) - 1|]
//     //let sC = Array.Parallel.mapi (fun i p ->
//     let sC = Array.mapi (fun i p ->
//                             let fp = finalDisparities.[i] |> int
//                             let mutable totalCost = 0.0f
//                             for (neighbourIdx, _) in p do
//                                 let fq = finalDisparities.[neighbourIdx] |> int
//                                 totalCost <- totalCost + smoothnessCosts.[fp, fq]
//                             totalCost
//                         ) messages |> Array.sum
//     dC + sC

let beliefpropagation parameters bpparameters =
    let dataCosts = Data.computeDataCosts parameters bpparameters.dataFunction
    let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters.smoothnessFunction
    let messages = initMessages parameters dataCosts
    let mutable indexIsEven = false
    for _i = 1 to bpparameters.iterations do
        updateMessages parameters.maximumDisparity dataCosts smoothnessCosts indexIsEven messages
        indexIsEven <- not indexIsEven

    // let findeps = computeFinalDisparities parameters dataCosts messages
    // findeps
    computeFinalDisparities parameters dataCosts messages

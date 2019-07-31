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

let inline truncateFH (minimum : ^a) (arr : ^a []) =
    Array.iteri (fun i x -> if minimum < x then arr.[i] <- minimum) arr

// This is pretty much directly lifted from F&H, both the paper and their sample code
let dtFH m =
    for fq = 1 to ((Array.length m) - 1) do
        m.[fq] <- min m.[fq] (m.[fq - 1] + Smoothness.C_FH)

    for fq = ((Array.length m) - 2) downto 0 do
        m.[fq] <- min m.[fq] (m.[fq + 1] + Smoothness.C_FH)

// This is intended to match eq. 2 in F & H 2006
//let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) (m1 : (int * float32 []) [] []) (m2 : (int * float32 []) [] []) =
let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) indexIsEven (m : (int * float32 []) [] []) =
// p and q are used below in accordance with Felzenswalb & Huttenlocher's notation
    // //Array.Parallel.iteri (fun i p -> // each pixel in the image
    // //Array.iteri (fun i p -> // each pixel in the image
    // let startIndex = if oddOrEven then 0 else 1
    // for i in startIndex..2..((Array.length m) - 1) do
    //     let p = m.[i]
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
    //) m

    // Using F&H's 'min convolution' style, where the message from p to q is computed as
    // min_{f_p} (V(f_p - f_q) + h(f_p)), where h(f_p) is
    // D_p(f_p) + the sum of the messages for f_p received from the other neighbours at t-1
    //Array.Parallel.iteri (fun i p -> // each pixel in the image
    //Array.iteri (fun i p -> // each pixel in the image
    // let startIndex = if oddOrEven then 0 else 1
    // let neighbourMessageSums = Array.zeroCreate (maxD + 1)
    // for i in startIndex..2..((Array.length m) - 1) do
    //     System.Array.Fill(neighbourMessageSums, 0.0f)
    //     let p = m.[i]
    //     for fp = 0 to maxD do
    //         for (_, (neighbourCosts : float32 [])) in p do
    //             neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]

    //     for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
    //         let hfp = Array.map2 (( - )) neighbourMessageSums neighbourCosts
    //         let minhfp = (Array.min hfp) + Smoothness.D_FH
    //         dtFH hfp
    //         truncateFH minhfp hfp
    //         let indexInNeighbour = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
    //         let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
    //         for fq = 0 to maxD do // each disparity label of q
    //             mneighbourcosts.[fq] <- min hfp.[fq] minhfp

    let startIndex = if indexIsEven then 0 else 1
    // let neighbourMessageSums = Array.zeroCreate (maxD + 1)
    let processPixel loopIndex _loop neighbourMessageSums =
        //fun loopIndex _loop neighbourMessageSums ->
        System.Array.Fill(neighbourMessageSums, 0.0f)
        let i =
            if indexIsEven then
                loopIndex * 2
            else
                loopIndex * 2 - 1
        let p = m.[i]
        for fp = 0 to maxD do
            for (_, (neighbourCosts : float32 [])) in p do
                neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]

        for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
            let hfp = Array.map2 (( - )) neighbourMessageSums neighbourCosts
            let minhfp = (Array.min hfp) + Smoothness.D_FH
            dtFH hfp
            truncateFH minhfp hfp
            let indexInNeighbour = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
            let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
            for fq = 0 to maxD do // each disparity label of q
                mneighbourcosts.[fq] <- min hfp.[fq] minhfp

        neighbourMessageSums

    let createThreadLocalArray = fun () -> Array.zeroCreate (maxD + 1)

    Threading.Tasks.Parallel.For(
        startIndex,
        ((Array.length m) / 2),
        createThreadLocalArray,
        processPixel,
        ignore
    ) |> ignore

    // for i in startIndex..2..((Array.length m) - 1) do
    //     System.Array.Fill(neighbourMessageSums, 0.0f)
    //     let p = m.[i]
    //     for fp = 0 to maxD do
    //         for (_, (neighbourCosts : float32 [])) in p do
    //             neighbourMessageSums.[fp] <- neighbourMessageSums.[fp] + neighbourCosts.[fp] + dataCosts.[i].[fp]

    //     for ((neighbourIdx : int), (neighbourCosts : float32 [])) in p do // each neighbour of the current pixel
    //         let hfp = Array.map2 (( - )) neighbourMessageSums neighbourCosts
    //         let minhfp = (Array.min hfp) + Smoothness.D_FH
    //         dtFH hfp
    //         truncateFH minhfp hfp
    //         let indexInNeighbour = Array.findIndex (fst >> ((=) i)) m.[neighbourIdx]
    //         let (_, mneighbourcosts) = m.[neighbourIdx].[indexInNeighbour]
    //         for fq = 0 to maxD do // each disparity label of q
    //             mneighbourcosts.[fq] <- min hfp.[fq] minhfp

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

let beliefpropagation parameters bpparameters =
    let dataCosts = Data.computeDataCosts parameters bpparameters.dataFunction
    let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters.smoothnessFunction
    let messages = initMessages parameters dataCosts
    let mutable indexIsEven = true // start with 0..2..etc
    for _i = 1 to bpparameters.iterations do
        updateMessages parameters.maximumDisparity dataCosts smoothnessCosts indexIsEven messages
        indexIsEven <- not indexIsEven

    computeFinalDisparities parameters dataCosts messages

// Currently just trying to focus on a basic version of min-sum 2D loopy belief propagation
// Using min-sum on the same basis that Felzenswalb & Huttenlocher substituted max-product with it

module BeliefPropagation

open Common
open Smoothness
open System

//[<Struct; IsByRefLike>]
type BPParameters<'a,'b> = {
    dataFunction: 'a -> 'a -> 'b
    dataCosts: 'b [] []
    smoothnessFunction: int -> int -> 'b
    smoothnessCosts: float32[,]
    iterations: int
}

//[<Struct; IsByRefLike>]
type Proxel<'a> = {
    index: int
    neighbours: int[]
    neighbourMessages: 'a[][]
    //dataCosts: 'a[]
}

// let computeDataCosts parameters bpparameters i =
//     let retArr = Array.zeroCreate (parameters.maximumDisparity + 1)
//     let leftIntensity = parameters.leftImage.[i]
//     for x = 0 to parameters.maximumDisparity do
//         retArr.[x] <- bpparameters.dataFunction leftIntensity (parameters.rightImage.[i - x])
//     retArr

// let computeDataCosts parameters bpparameters i =
//     Array.init (parameters.maximumDisparity + 1) (
//         fun j ->
//             bpparameters.dataFunction parameters.leftImage.[i] parameters.rightImage.[i - j]
//         )

// let computeNeighbours parameters i =
//     [|
//         (if i - 1 < 0 then None else Some(i - 1));
//         (if i + 1 >= parameters.width then None else Some(i + 1));
//         (if i - parameters.width < 0 then None else Some(i - parameters.width));
//         (if i + parameters.width >= parameters.height * parameters.width - 1 then None else Some(i + parameters.width));
//     |] |> Array.choose id

let computeNeighbours parameters i =
    [|
        if i - 1 < 0 then yield i - 1;
        if i + 1 >= parameters.width then yield i + 1;
        if i - parameters.width < 0 then yield i - parameters.width;
        if i + parameters.width >= parameters.height * parameters.width - 1 then yield i + parameters.width;
    |]

let makeProxels (parameters: Parameters<_>) bpparameters =
    // let bothImages = Array.zip parameters.leftImage parameters.rightImage
    // let dataCosts = Array.Parallel.map (fun x -> bpparameters.dataFunction <|| x) bothImages
    Array.Parallel.init (Array.length parameters.leftImage) (
        fun i ->
            let neigh = computeNeighbours parameters i
            {index = i
             //dataCosts = computeDataCosts parameters bpparameters i;
             neighbours = neigh;
             neighbourMessages = Array.zeroCreate (Array.length neigh)}
    )

let normalize messageArray =
    let divisor = Array.max messageArray
    Array.map (fun i -> i / divisor) messageArray


// Based on equation 2 on page 42 of Felzenswalb & Huttenlocher (2006)
let computeMessage parameters bpparameters neighbour proxel =
    let outgoingMessages = Array.zeroCreate parameters.maximumDisparity
    //let mutable minCost = Int32.MaxValue
    let mutable minCost = Single.MaxValue
    for i = 0 to parameters.maximumDisparity do
        //minCost <- Int32.MaxValue
        minCost <- Single.MaxValue
        for j = 0 to parameters.maximumDisparity do
            let messagesSum = (Array.sum proxel.neighbourMessages.[j]) - proxel.neighbourMessages.[j].[neighbour]
            //let smoothnessCost = bpparameters.smoothnessFunction i j
            let smoothnessCost = bpparameters.smoothnessCosts.[proxel.index, j]
            //let totalCost = messagesSum + smoothnessCost + proxel.dataCosts.[j]
            let totalCost = messagesSum + smoothnessCost + bpparameters.dataCosts.[proxel.index].[j]
            if totalCost < minCost then minCost <- totalCost
        outgoingMessages.[i] <- minCost

    normalize outgoingMessages

let computeMessages parameters bpparameters proxel =
    Array.mapi (computeMessage parameters bpparameters) proxel

let exchangeMessages proxels =
    ()

let computeFinalDisparities proxels =
    5.0f |> byte

// let beliefpropagation parameters bpparameters =
//     //let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters (Smoothness.pottsFloat32 Smoothness.LAMBDA_FH)
//     let dataCosts = Data.computeDataCosts parameters Data.absoluteDifference |> Array.map (Array.map float32)
//     let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters (Smoothness.potts Smoothness.LAMBDA_FH)
//     let proxels = makeProxels parameters {bpparameters with dataCosts = dataCosts; smoothnessCosts = smoothnessCosts}


//     for i = 1 to bpparameters.iterations do
//         exchangeMessages proxels // make all the proxels 'exchange messages' with one another

//     // Do the final computation to return the best-belief disparity
//     Array.map computeFinalDisparities proxels

let initMessages parameters =
    Array3D.zeroCreate (parameters.width * parameters.height) 4 parameters.maximumDisparity // 4 because I'm using the 4-neighbourhood

let updateMessages maxD dataCosts smoothnessCosts neighbours m1 m2 = // this function computes updates to the messages, using input data from m1 and storing it into m2, then swaps their pointers
    let findMin p q fq =
        let messageCosts = m1.[p, 0, fq] + m1.[p, 1, fq] + m1.[p, 2, fq] + m1.[p, 3, fq] - m1.[p, q, fq]
        let possibilities = Seq.init maxD
                                (fun fp ->
                                    smoothnessCosts.[fp, fq] + dataCosts.[p].[fp] + messageCosts
                                )
        Seq.min possibilities


    for p = 0 to (Array3D.length1 m1) do
        for q in neighbours.[p] do
            for fq = 0 to maxD do
                m2.[p,q,fq] <- findMin p q fq


let beliefpropagation parameters bpparameters =
    let neighbours = Array.init (parameters.width * parameters.height) (computeNeighbours parameters)
    let dataCosts = Data.computeDataCosts parameters Data.absoluteDifference // need to pass in the images, at least
    let smoothnessCosts = computeSmoothnessCosts parameters (Smoothness.potts Smoothness.LAMBDA_FH) // will need the maximum disparity for the smoothness costs
    let mutable messages1 = initMessages parameters // All messages will be 0 initially - need to work out how to represent the messages.  Big ol' array?
    let mutable messages2 = Array2D.copy messages1
    for i = 1 to bpparameters.iterations do
        updateMessages parameters.maximumDisparity messages1 messages2

    computeFinalDisparities
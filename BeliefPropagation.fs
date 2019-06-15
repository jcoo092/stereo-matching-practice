// Currently just trying to focus on a basic version of max-sum 2D loopy belief propagation

module BeliefPropagation

open Common
open System.Runtime.CompilerServices
open System

[<Struct; IsByRefLike>]
type BPParameters<'a,'b> = {
    dataFunction: 'a -> 'a -> 'b
    smoothnessFunction: int -> int -> 'b
    iterations: int
}

[<Struct; IsByRefLike>]
type Proxel<'a> = {
    neighbours: int[]
    neighbourMessages: 'a[][]
    dataCosts: 'a[]
}

// let computeDataCosts parameters bpparameters i =
//     let retArr = Array.zeroCreate (parameters.maximumDisparity + 1)
//     let leftIntensity = parameters.leftImage.[i]
//     for x = 0 to parameters.maximumDisparity do
//         retArr.[x] <- bpparameters.dataFunction leftIntensity (parameters.rightImage.[i - x])
//     retArr

let computeDataCosts parameters bpparameters i =
    Array.init (parameters.maximumDisparity + 1) (
        fun j ->
            bpparameters.dataFunction parameters.leftImage.[i] parameters.rightImage.[i - j]
        )

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
            {dataCosts = computeDataCosts parameters bpparameters i;
             neighbours = neigh;
             neighbourMessages = Array.zeroCreate (Array.length neigh)}
    )

let normalize messageArray =
    let divisor = Array.max messageArray
    Array.map (fun i -> i / divisor) messageArray


// Based on equation 2 on page 42 of Felzenswalb & Huttenlocher (2006)
let computeMessage parameters bpparameters neighbour proxel =
    let outgoingMessages = Array.zeroCreate parameters.maximumDisparity
    let mutable minCost = Int32.MaxValue
    for i = 0 to parameters.maximumDisparity do
        minCost <- Int32.MaxValue
        for j = 0 to parameters.maximumDisparity do
            let messagesSum = (Array.sum proxel.neighbourMessages.[j]) - proxel.neighbourMessages.[j].[neighbour]
            let smoothnessCost = bpparameters.smoothnessFunction i j
            let totalCost = messagesSum + smoothnessCost + proxel.dataCosts.[j]
            if totalCost < minCost then minCost <- totalCost
        outgoingMessages.[i] <- minCost

    normalize outgoingMessages

let computeMessages parameters bpparameters proxel =
    Array.mapi (computeMessage parameters bpparameters) proxel


let beliefpropagation parameters bpparameters =
    let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters (Smoothness.pottsFloat32 Smoothness.LAMBDA_FH)
    let proxels = makeProxels parameters bpparameters


    for i = 1 to bpparameters.iterations do
        exchangeMessages proxels // make all the proxels 'exchange messages' with one another

    // Do the final computation to return the best-belief disparity
    Array.map computeFinalDisparities proxels
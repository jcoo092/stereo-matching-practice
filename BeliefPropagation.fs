// Currently just trying to focus on a basic version of min-sum 2D loopy belief propagation
// Using min-sum on the same basis that Felzenswalb & Huttenlocher substituted max-product with it

module BeliefPropagation

open Common
open System
open MathNet.Numerics.LinearAlgebra

let standardNeighboursMatrix =
    matrix [|
        [|0.0f; 1.0f; 1.0f; 1.0f; 1.0f;|]; // Neighbour to the North
        [|1.0f; 0.0f; 1.0f; 1.0f; 1.0f;|]; // Neighbour to the East
        [|1.0f; 1.0f; 0.0f; 1.0f; 1.0f;|]; // Neighbour to the South
        [|1.0f; 1.0f; 1.0f; 0.0f; 1.0f;|]; // Neighbour to the West
    |]

[<Struct; Runtime.CompilerServices.IsByRefLike>]
type Proxel = {
    neighbourIndices: int []
    costMatrix: Matrix<float32>
    neighboursMatrix: Matrix<float32>
}

//[<Struct; Runtime.CompilerServices.IsByRefLike>]
type OptionProxel = {
    neighbourIndices: int ValueOption []
    costMatrix: Matrix<float32>  // Matrix that the costs for the messages are stored in
    neighboursMatrix: Matrix<float32>  // Matrix specifying the coefficients to apply to the various costs in the costMatrix
}

[<Struct; Runtime.CompilerServices.IsByRefLike>]
type BPParameters = {
    dataFunction: Byte -> Byte -> float32
    smoothnessFunction: int -> int -> float32
    iterations: int
}

let computeNeighbours parameters i =
    let x = i % parameters.width
    let y = i / parameters.width

    // let mutable x = 0
    // let y = Math.DivRem(i, parameters.width, &x)

    [|
        if y > 0 then yield i - parameters.width;
        if x > 0 then yield i - 1;
        if x < (parameters.width - 1) then yield i + 1;
        if y < (parameters.height - 1) then yield i + parameters.width;
    |]

let computeOptionNeighbours parameters i =
        let x = i % parameters.width
        let y = i / parameters.width

        // let mutable x = 0
        // let y = Math.DivRem(i, parameters.width, &x)

        let retArr = Array.create 4 ValueNone

        if y > 0 then retArr.[0] <- ValueSome(i - parameters.width)
        if x > 0 then retArr.[1] <- ValueSome(i - 1)
        if x < (parameters.width - 1) then retArr.[2] <- ValueSome(i + 1)
        if y < (parameters.height - 1) then retArr.[3] <- ValueSome(i + parameters.width)

        retArr

// let inline initMessages parameters =
//     //Array.Parallel.init (parameters.width * parameters.height) (
//     Array.init (parameters.width * parameters.height) (
//         fun i ->
//                 let neighbours = computeNeighbours parameters i
//                 Array.map (fun neighbour -> (neighbour, Array.zeroCreate (parameters.maximumDisparity + 1))) neighbours
//     )

let initOptionProxel parameters (dataCosts : float32 [][]) i =
    let neighbourOptionIndices = computeOptionNeighbours parameters i
    let neighboursMatrix =
        if Array.forall (ValueOption.isSome) neighbourOptionIndices then
            standardNeighboursMatrix
        else
            let nM = standardNeighboursMatrix.Clone()
            Array.iteri (
                fun i x ->
                    if ValueOption.isNone x then
                        nM.ClearColumn i
            ) neighbourOptionIndices
            nM

    let costMatrix = DenseMatrix.create 5 parameters.maximumDisparity Matrix.Zero
    costMatrix.SetRow(4, dataCosts.[i])

    {
        neighbourIndices = neighbourOptionIndices
        costMatrix = costMatrix
        neighboursMatrix = neighboursMatrix
    }

let inline initOptionProxels parameters dataCosts = Array.init (parameters.width * parameters.height - 1) (initOptionProxel parameters dataCosts)

let computeIndexInNeighbour =
    function
    | 0 -> 2 // Neighbour is to the North
    | 1 -> 3 // Neighbour is to the East
    | 2 -> 0 // Neighbour is to the South
    | 3 -> 1 // Neighbour is to the West
    | x -> failwith (sprintf "Invalid neighbourIndexArrayPosition: %d" x)

let inline normalizeCostArray arr =
    let mini = Array.min arr
    Array.iteri (fun i value -> arr.[i] <- value / mini) arr

let inline computeAndSendNewMessages parameters (smoothnessCosts : float32 [,]) (proxels : OptionProxel[]) proxelIndex proxel =
    let outgoingMessages = proxel.neighboursMatrix.Multiply(proxel.costMatrix)

    // Array.iteri (
    //     fun i v ->
    //     match v with
    //     | ValueSome(neighbourIndex) ->
    //         let costsWithoutSmoothness = outgoingMessages.AsRowArrays().[neighbourIndex]
    //         printfn "got past the declaration of costwithoutsmoothness!"
    //         let scratchSpaceArray = Array.zeroCreate parameters.maximumDisparity
    //         let finalCosts = Array.init parameters.maximumDisparity
    //                             (fun j ->
    //                             Array.iteri (
    //                                 fun k _ ->
    //                                 scratchSpaceArray.[k] <- costsWithoutSmoothness.[k] + smoothnessCosts.[j,k]
    //                             ) scratchSpaceArray
    //                             Array.min scratchSpaceArray
    //                             )
    //         normalizeCostArray finalCosts
    //         proxels.[neighbourIndex].costMatrix.SetRow(computeIndexInNeighbour i, finalCosts)
    //     | ValueNone -> ()
    // ) proxel.neighbourIndices

    Seq.iteri (
        fun i (costsWithoutSmoothness : Vector<float32>) ->
            match proxel.neighbourIndices.[i] with
            | ValueSome(neighbourIndex) ->
                let scratchSpaceArray = Array.zeroCreate parameters.maximumDisparity
                let finalCosts = Array.init parameters.maximumDisparity
                                    (fun j ->
                                    Array.iteri (
                                        fun k _ ->
                                        scratchSpaceArray.[k] <- costsWithoutSmoothness.[k] + smoothnessCosts.[j,k]
                                    ) scratchSpaceArray
                                    Array.min scratchSpaceArray
                                    )
                normalizeCostArray finalCosts
                //printfn "finished normalizing costs!"
                proxels.[neighbourIndex].costMatrix.SetRow(computeIndexInNeighbour i, finalCosts)
            | ValueNone -> ()
    ) (outgoingMessages.EnumerateRows())
    // printfn "Just finished processing proxel number %d" proxelIndex

// // This is intended to match eq. 2 in F & H 2006
// let updateMessages maxD (dataCosts : float32 [][]) (smoothnessCosts : float32 [,]) (m1 : (int * float32 []) [] []) (m2 : (int * float32 []) [] []) =
// // this function computes updates to the messages using data in m1, and stores it back to m2
// // p and q are used below in accordance with Felzenswalb & Huttenlocher's notation
//     //Array.Parallel.iteri (fun i p -> // each pixel in the image
//     Array.iteri (fun i p -> // each pixel in the image
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
//     //Array.Parallel.iter (fun p ->
//     Array.iter (fun p ->
//                         Array.iter (fun (_, neighbourCosts) -> normalizeCostArray neighbourCosts) p
//     ) m2

// let computeFinalDisparities parameters (dataCosts : float32 [][]) (messages : (int * float32 []) [] []) =
//     //Array.Parallel.mapi (fun i p ->
//     Array.mapi (fun i p ->
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

let computeFinalDisparities (proxels : OptionProxel[]) =
    Array.map (fun p -> Vector.minIndex (p.costMatrix.ColumnSums()) |> byte) proxels

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

let getOddOrEvenProxels proxels oddOrEven =
    seq {
        let startingIndex =
            if oddOrEven then
                1
            else
                0
        for i in startingIndex..2..(Array.length proxels - 1) do
            yield proxels.[i]
    } |> Seq.skip 50000

let beliefpropagation parameters bpparameters =
    let dataCosts = Data.computeDataCosts parameters bpparameters.dataFunction
    let smoothnessCosts = Smoothness.computeSmoothnessCosts parameters bpparameters.smoothnessFunction
    //let mutable messages1 = initMessages parameters // All messages will be 0 initially
    //let mutable messages2 = initMessages parameters
    let proxels = initOptionProxels parameters dataCosts
    let mutable oddOrEven = false
    for _i = 1 to bpparameters.iterations do
        //updateMessages parameters.maximumDisparity dataCosts smoothnessCosts messages1 messages2
        Seq.iteri (computeAndSendNewMessages parameters smoothnessCosts proxels) (getOddOrEvenProxels proxels oddOrEven)
        // let temp = messages1
        // messages1 <- messages2
        // messages2 <- temp
        oddOrEven <- not oddOrEven

    //let findeps = computeFinalDisparities parameters dataCosts messages1
    let findeps = computeFinalDisparities proxels
    //let finenergy = computeEnergy dataCosts smoothnessCosts messages1 findeps
    // printfn "Final energy is: %f" (finenergy / float32 parameters.totalPixels)
    findeps

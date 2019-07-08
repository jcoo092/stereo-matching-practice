// Currently just trying to focus on a basic version of min-sum 2D loopy belief propagation
// Using min-sum on the same basis that Felzenswalb & Huttenlocher substituted max-product with it

module BeliefPropagation

open Common
open System

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

let initMessages parameters (dataCosts : float32[][]) =
    Array.Parallel.init (parameters.width * parameters.height) (
        fun i ->
            let message = Array2D.create 5 parameters.maximumDisparity 0.0f
            for j = 0 to (parameters.maximumDisparity - 1) do
                message.[int Direction.Data, j] <- dataCosts.[i].[j]
            {msg = message; bestDisparity = 0}
    )

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
    let newMsg = Array2D.create 1 parameters.maximumDisparity 0.0f
    let width = parameters.width
    for i = 0 to (parameters.maximumDisparity - 1) do
        let mutable minVal = Single.MaxValue
        for j = 0 to (parameters.maximumDisparity - 1) do
            let mutable p = smoothnessCosts.[i, j] + messages.[x + y * width].msg.[int Direction.Data, j]

            if direction <> Direction.West then
                p <- p + messages.[x + y * width].msg.[int Direction.West, j]
            if direction <> Direction.East then
                p <- p + messages.[x + y * width].msg.[int Direction.East, j]
            if direction <> Direction.North then
                p <- p + messages.[x + y * width].msg.[int Direction.North, j]
            if direction <> Direction.South then
                p <- p + messages.[x + y * width].msg.[int Direction.South, j]

            minVal <- min minVal p

        newMsg.[0,i] <- minVal

    let minVal =
        let mini = Array.min newMsg.[0,*]
        if mini < 1.0f then
            1.0f
        else
            mini
    Array2D.iteri (fun i j x -> newMsg.[i,j] <- x / minVal) newMsg

    match direction with
    | Direction.West ->
        Array2D.blit newMsg 0 0 messages.[y*width + (x - 1)].msg (int Direction.East) 0 1 parameters.maximumDisparity
    | Direction.East ->
        Array2D.blit newMsg 0 0 messages.[y*width + (x + 1)].msg (int Direction.West) 0 1 parameters.maximumDisparity
    | Direction.North ->
        Array2D.blit newMsg 0 0 messages.[(y - 1)*width + x].msg (int Direction.South) 0 1 parameters.maximumDisparity
    | Direction.South ->
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
        bp parameters smoothnessCosts messages Direction.East
        bp parameters smoothnessCosts messages Direction.West
        bp parameters smoothnessCosts messages Direction.North
        bp parameters smoothnessCosts messages Direction.South

        let energy = MAP parameters smoothnessCosts messages

        printfn "iteration %d, energy = %f" i energy

    Array.map (fun m -> m.bestDisparity |> byte) messages
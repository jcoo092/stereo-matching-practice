#load "Common.fs"
#load "Data.fs"
#load "Smoothness.fs"
#load "BeliefPropagation.fs"
#time

let leftInput = [|
    127uy; 126uy; 123uy; 129uy;
    55uy; 49uy; 53uy; 58uy;
    48uy; 47uy; 50uy; 52uy;
|]

let rightInput = [|
    127uy; 127uy; 124uy; 130uy;
    53uy; 49uy; 54uy; 57uy;
    48uy; 201uy; 203uy; 199uy;
|]

//printfn "\nDefined inputs\n"

let matchingParameters : Common.Parameters = {
        leftImage = leftInput
        rightImage = rightInput
        width = 4
        height = 3
        totalPixels = 12
        windowEdgeSize = 3
        maximumDisparity = 4
        zeroMean = true
}

let bpparameters : BeliefPropagation.BPParameters = {
    dataFunction = (Data.FHTruncatedLinear Smoothness.LAMBDA_FH Smoothness.TAU_FH)
    smoothnessFunction = (Smoothness.truncatedLinear Smoothness.D_FH)
    iterations = 4
}

//printfn "\nMade it past initial declarations\n"

let dataCosts = Data.computeDataCosts matchingParameters bpparameters.dataFunction
//printfn "%A" dataCosts
let smoothnessCosts = Smoothness.computeSmoothnessCosts matchingParameters bpparameters.smoothnessFunction
//printfn "smooth dims are %d & %d" (Array2D.length1 smoothnessCosts) (Array2D.length2 smoothnessCosts)
let mutable messages1 = BeliefPropagation.initMessages matchingParameters // All messages will be 0 initially
//printfn "%A" messages1
let mutable messages2 = BeliefPropagation.initMessages matchingParameters
for _i = 1 to bpparameters.iterations do
        BeliefPropagation.updateMessages matchingParameters.maximumDisparity dataCosts smoothnessCosts messages1 messages2
        let temp = messages1
        messages1 <- messages2
        messages2 <- temp

printfn "m1:\n%A\n" messages1
printfn "m2:\n%A\n" messages2
printfn "%A" (messages1 = messages2)
module Smoothness

// d here matches to 'd' used in Felzenswalb & Huttenlocher (2006) - it is a truncation parameter
let inline truncatedLinear d a b =
    (min (Data.absoluteDifference a b |> float32) d)

let inline computeSmoothnessCosts (parameters : Common.Parameters) (smoothnessFunction : int -> int -> single) =
    Array2D.init (parameters.maximumDisparity + 1) (parameters.maximumDisparity + 1) smoothnessFunction
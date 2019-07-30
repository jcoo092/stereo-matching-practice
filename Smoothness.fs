module Smoothness

[<Literal>]
let LAMBDA_FH = 0.07f
[<Literal>]
let TAU_FH = 15.0f
[<Literal>]
let C_FH = 1.0f
[<Literal>]
let D_FH = 1.7f
[<Literal>]
let EPSILON_JAMES = 7.62939453125e-5f // Very small value, equivalent to 1e-17 (if I remember rightly)

let inline potts lambda a b =
    if a = b then
        LanguagePrimitives.GenericZero
    else
        lambda

let inline pottsFloat32 (lambda: float32) (a : float32) b =
    if System.Math.Abs(a - b) < EPSILON_JAMES then
        0.0f
    else
        lambda

// d here matches to 'd' used in Felzenswalb & Huttenlocher (2006) - it is a truncation parameter
let inline truncatedLinear d a b =
    (min (Data.absoluteDifference a b |> float32) d)

let inline truncatedQuadratic lambda d a b =
    (min (Data.squaredDifference a b) d) |> float32 |> (*) lambda

let inline computeSmoothnessCosts (parameters : Common.Parameters) (smoothnessFunction : int -> int -> single) =
    Array2D.init (parameters.maximumDisparity + 1) (parameters.maximumDisparity + 1) smoothnessFunction
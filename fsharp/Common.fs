module Common

open System

type LeftOrRight = | Left | Right

[<Struct;>]
type Parameters = {
    leftImage: byte []
    rightImage: byte []
    width: int
    height: int
    totalPixels: int
    windowEdgeSize: int
    maximumDisparity: int
    zeroMean: bool
}

[<Literal>]
let LAMBDA_FH = 0.07f
[<Literal>]
let TAU_FH = 15.0f
[<Literal>]
let C_FH = 1.0f
[<Literal>]
let D_FH = 1.7f
[<Literal>]
let CUSTOM_EPSILON = 7.62939453125e-5f // Very small value, equivalent to 1e-17 (if I remember rightly)

// Taken from 'Average of Integers', listed at http://aggregate.org/MAGIC/#Average%20of%20Integers accessed on 14 June 2019
let inline twoParameterMean a b =
    (a &&& b) + ((a ^^^ b) >>> 1)

let inline saturatingSubtraction minuend subtrahend =
    if subtrahend > minuend then
        LanguagePrimitives.GenericZero
    else
        minuend - subtrahend

let float32Equality a b =
    let difference =
        if a < b then
            b - a
        else
            a - b
    difference < Single.Epsilon

let buildArraySlices parameters usedImage =
    let subjectArray = match usedImage with
                        | Left -> parameters.leftImage
                        | Right -> parameters.rightImage
    Array.chunkBySize parameters.width subjectArray

let argminFloat32Array arr =
        let min = Array.min arr
        Array.findIndex (float32Equality min) arr
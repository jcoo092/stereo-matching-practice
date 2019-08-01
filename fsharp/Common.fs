module Common

open System
//open FSharp.Span.Utils

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

let buildArraySegments parameters usedImage =
    let subjectArray = match usedImage with
                        | Left -> parameters.leftImage
                        | Right -> parameters.rightImage
    let width = parameters.width
    let slices = Array.zeroCreate parameters.height
    for i in 0..(parameters.height - 1) do
        slices.[i] <- ArraySegment(subjectArray, i * width, width)
    slices

let argminArray arr =
    let min = Array.min arr
    Array.findIndex ((=) min) arr

let argminFloat32Array arr =
        let min = Array.min arr
        //printfn "The current min is: %f" min
        Array.findIndex (float32Equality min) arr

// The below is by Tomas Petricek
// Directly copied from http://fssnip.net/1R/title/Take-every-Nth-element-of-sequence
// [snippet:Efficient version using enumerator]
let everyNth n (input:seq<_>) =
  seq { use en = input.GetEnumerator()
        // Call MoveNext at most 'n' times (or return false earlier)
        let rec nextN n =
          if n = 0 then true
          else en.MoveNext() && (nextN (n - 1))
        // While we can move n elements forward...
        while nextN n do
          // Retrun each nth element
          yield en.Current }
// [/snippet]
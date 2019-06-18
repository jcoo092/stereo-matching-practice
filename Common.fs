module Common

open System
open FSharp.Span.Utils

type LeftOrRight = | Left | Right

//[<Struct; IsByRefLike>]
type Parameters<'a> = {
    //leftImage: byte []
    //leftImage: ReadOnlySpan<byte>
    leftImage: 'a []
    //rightImage: byte []
    //rightImage: ReadOnlySpan<byte>
    rightImage: 'a []
    width: int
    height: int
    windowEdgeSize: int
    maximumDisparity: int
    zeroMean: bool
}

// Taken from 'Average of Integers', listed at http://aggregate.org/MAGIC/#Average%20of%20Integers accessed on 14 June 2019
let inline twoParameterMean a b =
    (a &&& b) + ((a ^^^ b) >>> 1)

let inline saturingSubtraction minuend subtrahend =
    if subtrahend > minuend then
        LanguagePrimitives.GenericZero
    else
        minuend - subtrahend

let buildSlices parameters usedImage = // leftImage is a boolean, specifying whether the left (true) or right (false) image should be used
    //let slicer = ReadOnlyMemory(if useLeftImage then parameters.leftImage else parameters.rightImage)
    //let slicer = if usedImage then parameters.leftImage else parameters.rightImage
    let subjectArray = match usedImage with
                        | Left -> parameters.leftImage
                        | Right -> parameters.rightImage
                        |> ReadOnlySpan.ofArray
    let width = parameters.width
    let slices = Array.zeroCreate parameters.height
    for i in 0..(parameters.height - 1) do
        //slices.[i] <- subjectArray.Slice(i * width, width)
        slices.[i] <- ReadOnlySpan.slice (i * width) width subjectArray
    slices

let buildArraySlices parameters usedImage =
    let subjectArray = match usedImage with
                        | Left -> parameters.leftImage
                        | Right -> parameters.rightImage
    Array.chunkBySize parameters.width subjectArray

let buildArraySegments parameters usedImage =
    //let subjectArray = if useLeftImage then parameters.leftImage else parameters.rightImage
    let subjectArray = match usedImage with
                        | Left -> parameters.leftImage
                        | Right -> parameters.rightImage
    let width = parameters.width
    let slices = Array.zeroCreate parameters.height
    for i in 0..(parameters.height - 1) do
        slices.[i] <- ArraySegment(subjectArray, i * width, width)
    slices
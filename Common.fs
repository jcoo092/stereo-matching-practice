module Common

open System
open System.Runtime.CompilerServices
open FSharp.Span.Utils
open System.Runtime.CompilerServices

type LeftOrRight = | Left | Right

[<Struct; IsByRefLike>]
type Parameters = {
    leftImage: byte []
    //leftImage: ReadOnlySpan<byte>
    rightImage: byte []
    //rightImage: ReadOnlySpan<byte>
    width: int
    height: int
    windowEdgeSize: int
    maximumDisparity: int
    zeroMean: bool
}

let inline squaredDifference a b = pown (a - b) 2
let inline absoluteDifference a b = abs (a - b)

let inline arraysSquaredDifference (a: ^a[]) (b: ^a[]) i d =
    let a' = a.[i]
    let b' = b.[i - d]
    squaredDifference a' b'

let inline arraysAbsoluteDifference (a: ^a[]) (b: ^a[]) i d =
    let a' = a.[i]
    let b' = b.[i - d]
    absoluteDifference a' b'

let inline slicesSquaredDifference (a: ReadOnlyMemory< ^a >) (b: ReadOnlyMemory< ^a >) i d =
    let a' = a.Span.[i]
    let b' = b.Span.[i - d]
    squaredDifference a' b'

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
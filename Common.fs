module Common

open System
open System.Runtime.CompilerServices
open FSharp.Span.Utils

type LeftOrRight = | Left | Right

[<Struct; IsByRefLike>]
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

let inline squaredDifference a b = pown (a - b) 2
let inline absoluteDifference a b = abs (a - b)

let inline arraysSquaredDifference (a: ^a []) (b: ^a []) i d =
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

// This function is directly lifted from A Pixel Dissimilarity Measure That Is Insensitive to Image Sampling (1998) by Birchfield & Tomasi
// Interestingly, in their notation they seem to operate on the basis of a functional representation of the image
// That is, they have a function that takes a given coordinate, and returns the corresponding intensity
// For the below, ln is the I^-_L and lp is the I^+_L, while l is I_L(x_L), and much the same for the right image values
let inline birchfieldTomasi ln l lp rn r rp = 
    // For the below, a corresponds to I_L(x_l), b to I_R(x_r), c to x_{r - 1}, d to x_{r + 1}
    let inline computeDBar a b c d = 
        let irminus = twoParameterMean b c
        let irplus = twoParameterMean b d
        let imin = min b (min irminus irplus)
        let imax = max b (max irminus irplus)
        max LanguagePrimitives.GenericZero (max (saturingSubtraction a imax) (saturingSubtraction imin a))

    let dbarleft = computeDBar l r rn rp
    let dbarright = computeDBar r l ln lp
    min dbarleft dbarright

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
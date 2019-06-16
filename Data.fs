module Data

open System
open Common

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

let computeDataCosts parameters dataCostFunction =
    let data = Array.zeroCreate (parameters.width * parameters.height)
    for x = 0 to parameters.width do
        for y = 0 to parameters.height do
            let leftIdx = x + y * parameters.height
            let lowerBound = System.Math.Clamp(leftIdx - parameters.maximumDisparity, 0, leftIdx - parameters.maximumDisparity)
            let currentPixelData = Array.zeroCreate (leftIdx - lowerBound)
            for d = parameters.maximumDisparity downto lowerBound do
                currentPixelData.[d] <- dataCostFunction parameters.leftImage.[leftIdx] parameters.rightImage.[leftIdx - d]
            data.[leftIdx] <- currentPixelData
    data
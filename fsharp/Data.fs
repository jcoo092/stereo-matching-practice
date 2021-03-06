module Data

open System
open Common

let inline squaredDifference a b = pown (a - b) 2
let inline absoluteDifference a b = abs (a - b)

let inline carefulAbsoluteDifference a b =
    let bigger = max a b
    let smaller = min a b
    bigger - smaller

let inline manualAbsoluteDifference a b =
    let retVal = if a < b then b - a else a - b
    float32 retVal

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
        max LanguagePrimitives.GenericZero (max (saturatingSubtraction a imax) (saturatingSubtraction imin a))

    let dbarleft = computeDBar l r rn rp
    let dbarright = computeDBar r l ln lp
    min dbarleft dbarright

// TODO:  Finish this!
let inline btDifference a b = 5

let inline FHTruncatedLinear lambda tau a b =
    lambda * (min (manualAbsoluteDifference a b) tau)

let computeDataCosts (parameters: Common.Parameters) (dataCostFunction: byte -> byte -> single) =
    let data =
        Array.zeroCreate (parameters.width * parameters.height)

    for x = 0 to (parameters.width - 1) do
        for y = 0 to (parameters.height - 1) do
            let leftIdx = x + y * parameters.width

            let currentPixelData =
                Array.zeroCreate (parameters.maximumDisparity + 1)

            for d = parameters.maximumDisparity downto 0 do
                currentPixelData.[d] <-
                    if x - d < 0 then
                        dataCostFunction 255uy 0uy
                    else
                        dataCostFunction parameters.leftImage.[leftIdx] parameters.rightImage.[leftIdx - d]

            data.[leftIdx] <- currentPixelData

    data

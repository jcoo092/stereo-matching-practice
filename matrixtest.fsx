#I @"C:\Users\jcoo092\.nuget\packages\mathnet.numerics\4.8.1\lib\netstandard2.0\"
#I @"C:\Users\jcoo092\.nuget\packages\mathnet.numerics.fsharp\4.8.1\lib\netstandard2.0\"

#r "MathNet.Numerics"
#r "MathNet.Numerics.FSharp"

open MathNet.Numerics.LinearAlgebra

let productMatrix = matrix [|[|0.0f; 1.0f; 1.0f; 1.0f; 1.0f|];
                            [|1.0f; 0.0f; 1.0f; 1.0f; 1.0f;|];
                            [|1.0f; 1.0f; 0.0f; 1.0f; 1.0f;|];
                            [|1.0f; 1.0f; 1.0f; 0.0f; 1.0f|]|]


let valuesMatrix : float32 Matrix = DenseMatrix.random 5 16 (MathNet.Numerics.Distributions.Normal.WithMeanStdDev(10.0, 3.0))

//valuesMatrix.[4, 0..15] <- vector [|5.0f; 5.0f; 5.0f; 5.0f; 5.0f|]
valuesMatrix.[4, 0..15] <- DenseVector.init 16 (fun i -> 20.0f + (float32 i))

let finalMatrix = productMatrix * valuesMatrix

printfn "%A" finalMatrix

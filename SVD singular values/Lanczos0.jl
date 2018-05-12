
try
    include("src/LogGenerator.jl")
    include("src/Matrices.jl")
    using GenericSVD
catch ex
    print("Что-то пошло не так: ", ex)
end

len = length(alpha)

    k = div(len, 2)


            alpha[k] = alpha[k] - abs(betta[k])
            alpha[k+1] = alpha[k+1] - abs(betta[k])


        matr = toDense(alpha[1:k], betta[1:k-1])
        X1 = det(x*eye(k) - matr)

        matr = toDense(alpha[k+1:end], betta[k+1:end])
        X2 = det(x*eye(len-k) - matr)

        matr = toDense(alpha[1:k-1], betta[1:k-2])
        X3 = det(x*eye(k-1) - matr)

        matr = toDense(alpha[k+2:end], betta[k+2:end])
        X4 = det(x*eye(len-k-1) - matr)

    X0 = X1*X2 + abs(betta[k])*X1*X4 + abs(betta[k])*X2*X3

function changealpha(alpha)
    alpha[1] = 4
end

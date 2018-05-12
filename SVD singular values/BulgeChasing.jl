include("Lanczos.jl")

function givensMatrix(size, a, b, numDiagonalElement=1, blockSize=1)

    G = eye(size)
    r = sqrt(a^2+b^2)
    c = a/r
    s = -b/r
    G[numDiagonalElement, numDiagonalElement] = c
    G[numDiagonalElement, numDiagonalElement+blockSize] = -s
    G[numDiagonalElement+blockSize, numDiagonalElement] = s
    G[numDiagonalElement+blockSize, numDiagonalElement+blockSize] = c

    return G
end

function givensMatr(size, a, b, row=1, column=1, blockSize=1)

    G = eye(size)
    r = sqrt(a^2+b^2)
    c = a/r
    s = -b/r
    G[row, column] = c
    G[row, column+blockSize] = -s
    G[row+blockSize, column] = s
    G[row+blockSize, column+blockSize] = c

    return G
end



function applyGivensOnMatrix(A)

    n = size(A)[1]

    G = givensMatrix(n, A[1,1], A[2,1])
    A[1:3,1:3] = G * A[1:3,1:3] * G'

    for i = 1:1:(n-2)
        print("\n",size(A))
        G = givensMatrix(n, A[i+1,i], A[i+2,i], 2, 3)
        A[i:i+2,i:i+2] = G * A[i:i+2,i:i+2] * G'
    end

    return A
end

function bulgeChasing(A)

    n = size(A)[1]

    k=1
    while (A[k+1,k]<1e-10)
        k = k+1
        #(k>n-1)?continue:0
    end
    print("valK:",k)
    G = givensMatrix(n, A[k,k], A[k+1,k], k)
    A = G*A*G'
    print("ASDF ",A[k,k])
    #=for i in k:1:n-2

        G = givensMatrix(n, A[i+1,i], A[i+2,i], i+1)
        A = G*A*G'

    end=#

    return [A,G]
end


#alphas = [1 3 0.000007 9999 5 8]
#bettas = [2 500 9 1 7]
Af = toDense(alphas, bettas)
#=
G = givensMatrix(6, Af[1,1], Af[2,1], 1, 1)
A = G * Af * G'

G2 = givensMatrix(6, A[1,1], A[3,1], 1, 2)
A2 = G2 * A * G2'

G3 = givensMatrix(6, A2[1,1], A2[2,1], 1, 1)
A3 = G3 * A2 * G3'=#

#Af = A2



G1 = givensMatr(6, Af[1,1], Af[2,1], 1, 1,1)
A1 = G1*Af*G1'


G2 = givensMatr(6, A1[2,1], A1[3,1], 2, 2, 1)
A2 = G2 * A1 * G2'

G3 = givensMatrix(6, A2[3,2], A2[4,2], 3)

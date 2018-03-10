function copyDiagonal(A)
    n = size(A)[1]
    for k in 2:1:n
        A[k-1,k] = A[k,k-1]
    end
    return A
end

""" Algorithm from Lioyd - Numerical Line page 277 (for Symmetric) """
function Lanczos(A)

    n = size(A)[1]

    T = zeros(n+1,n)
    b = rand(n)
    q = b/norm(b)
    q_prev = 0
    a = 0
    b = 0
    for k in 1:1:n
        v = A*q
        T[k,k] = q'*v

        v = v - b*q_prev - T[k,k]*q
        b = norm(v)
        T[k+1,k] = b
        q_prev = q
        q = v/b
    end

    return copyDiagonal(T[1:n,:])
end



print("test Lanczos\n")

A = rand(100,30)
symMatrix = A'*A
res = Lanczos(symMatrix)
nev = eigvals(res)-eigvals(symMatrix)

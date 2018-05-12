
try
    include("src/LogGenerator.jl")
    include("src/Matrices.jl")
    using GenericSVD
catch ex
    print("Что-то пошло не так: ", ex)
end

""" Algorithm from LioydTtrefethen - Numerical Line page 277 (Modified for not Symmetric) """
function lanczos(A)

    m = size(A)[1]
	n = size(A)[2]
    ### alpha - diagonal elements
   	alpha = Array{BigFloat}(n)

   	### betta - subdiagonal elements
   	betta = Array{BigFloat}(n-1)

    b = rand(n)
    q = b/norm(b)
    q_prev = 0
    a = 0
    b = 0
    for k in 1:1:n-1
        v = A'*(A*q)
        alpha[k] = q'*v

        v = v - b*q_prev - alpha[k]*q
        b = norm(v)
        betta[k] = b
        q_prev = q
        q = v/b
    end

    v = A'*(A*q)
    alpha[n] = q'*v

    return [alpha, betta]
end

""" Алгоритм Ланцоша возвращающий плотную трехдиагональную матрицу """
function lanczosFull(A)

    m = size(A)[1]
	n = size(A)[2]

    Ar = copy(A)

    b = rand(n)
    q = b/norm(b)
    q_prev = 0
    a = 0
    b = 0
    for k in 1:1:n-1
        v = A'*(A*q)
        Ar[k,k] = q'*v

        v = v - b*q_prev - Ar[k,k]*q
        b = norm(v)
        Ar[k+1,k] = b
        Ar[k,k+1] = b
        q_prev = q
        q = v/b
    end

    v = A'*(A*q)
    Ar[n,n] = q'*v

    return Ar
end

function toDense(alpha, betta)

    len = length(alpha)
    A = zeros(len, len)

    for i in 1:1:(len-1)
        A[i,i]=alpha[i]
        A[i+1,i] = A[i,i+1] = betta[i]
    end
    A[end,end] = alpha[end]

    return A
end

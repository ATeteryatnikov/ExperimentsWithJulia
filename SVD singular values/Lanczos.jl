
try
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

    b = big.(ones(n))
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

function toDense(alpha, betta)

    len = length(alpha)
    A = fill(big(0e0), len, len)

    for i in 1:1:(len-1)
        A[i,i]=alpha[i]
        A[i+1,i] = A[i,i+1] = betta[i]
    end
    A[end,end] = alpha[end]

    return A
end

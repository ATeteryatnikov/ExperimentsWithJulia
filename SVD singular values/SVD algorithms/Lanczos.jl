
""" Algorithm from Lioyd - Numerical Line page 277 (Modified for not Symmetric) """
function Lanczos(A)

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

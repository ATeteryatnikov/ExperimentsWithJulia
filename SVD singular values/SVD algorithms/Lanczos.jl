
function getGodunovMatrix(godunov_matrix_dim)

n = godunov_matrix_dim
A = fill(big(0e0),n,n)  # Матрица A размерности n*n типа BigFloat задается
                       # функцией fill() и заполняется нулевыми элементами.
B = Array(BigFloat, n) # Вектор B – массив размерности n типа BigFloat.

for i in 1:1:n
  A[i,i] = 7/5  # При заполнении элементов массива, необходимо сначала
                   # преобразовать целое число в тип BigFloat, и только
                   # после этого, делить его на другое число.
end

for i in 1:1:(n-1)
 A[i,i+1] = 11/3
end

for i in 1:1:(n-2)   #Заполнение вектора B
 B[i] = (152*i+118)
 B[i] = B[i]/(15*(2*i+1)*(2*i+3))
end
B[n] = 7
B[n] = B[n]/(5*(2*n+1))
B[n-1] = (152*n-34)
B[n-1] = B[n-1]/(15*(2*n-1)*(2*n+1))

real_result = Array(BigFloat, n)
for i in 1:1:(n-2)
  real_result[i] = 1/(2*i+1)
end
real_result[n-1] = 1/(2*n-1)
real_result[n] = 1/(2*n+1)

return [A, B, real_result]
end



""" Algorithm from Lioyd - Numerical Line page 277 (Modified for not Symmetric) """
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

function givensMatrix(a, b, len, a_i=1, b_i=2)
    G = eye(len)
    r = sqrt(a^2+b^2)
    c = a/r
    s = -b/r
    G[a_i, a_i] = c
    G[a_i, b_i] = -s
    G[b_i, a_i] = s
    G[b_i, b_i] = c
    return G
end


function applyGivensOnMatrix(A)

    n = size(A)[1]
    G2 = eye(n)
    for i = 1:1:(n-1)
        G = givensMatrix(A[i,i], A[i+1,i], n, i, i+1)
        A = G * A
        G2 = G*G2
    end

    return [A,G2]
end

""" Вычисление полинома в точке x по схеме 4.6 """
function FDACA(x, alpha, betta, lenBlock, s="", mode="T")

    print("Level: ", s, " Mode : ", mode, "\n")
    level = length(s)
    len = length(alpha)
    print(len)
    #print("WRONG:  ", (len<=lenBlock))
    if (len <= lenBlock)
        #print("LENGTH RIGHT: ", len,"!\n")
        matr = toDense(alpha, betta)
        return det(matr - x*eye(len))
    else
        k = div(len, 2)
        #print("\n",k, "is K \n  ")

        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k+1] = alpha[k+1] - abs(betta[k])

        X1 = X2 = X3 = X4 = 1

        if (mode == "T")
            print("call : T")
            if (k-1 >= 0)
                X1 = FDACA(x, alpha[1:k], betta[1:k-1], lenBlock, string(s,"0"), "T")
                #matr = toDense(alpha[1:k], betta[1:k-1])
                #X1 = det(matr - x*eye(k))
                #print("len alpha = ", length(alpha[1:k]), "\n")
            end
            if (len-(k+1) >= 0)
                X2 = FDACA(x, alpha[k+1:end], betta[k+1:end], lenBlock, string(s,"1"), "T")
                #matr = toDense(alpha[k+1:end], betta[k+1:end])
                #X2 = det(matr - x*eye(len-k))
                #print("len alpha = ", length(alpha[k+1:end]), "\n")
            end
            if (k-2 >= 0)
                X3 = FDACA(x, alpha[1:k-1], betta[1:k-2], lenBlock, string(s,"0"), "S")
                #matr = toDense(alpha[1:k-1], betta[1:k-2])
                #X3 = det(matr - x*eye(k-1))
                #print("len alpha = ", length(alpha[1:k-1]), "\n")
            end
            if (len-(k+2) >= 0)
                X4 = FDACA(x, alpha[k+2:end], betta[k+2:end], lenBlock, string(s,"1"), "R")
                #matr = toDense(alpha[k+2:end], betta[k+2:end])
                #X4 = det(matr - x*eye(len-k-1))
                #print("len alpha = ", length(alpha[k+2:end]), "\n")
            end
            print("call : T", s, " end \n")
        end

        if (mode == "S")
            print("call : S")
            if (k-1 >= 0)
                X1 = FDACA(x, alpha[1:k], betta[1:k-1], lenBlock, string(s,"0"), "T")
            end
            if ((len-1)-(k+1) >= 0)
                X2 = FDACA(x, alpha[k+1:end-1], betta[k+1:end-1], lenBlock, string(s,"1"), "S")
            end
            if (k-2 >= 0)
                X3 = FDACA(x, alpha[1:k-1], betta[1:k-2], lenBlock, string(s,"0"), "S")
            end
            if ((len-1)-(k+2) >= 0)
                X4 = FDACA(x, alpha[k+2:end-1], betta[k+2:end-1], lenBlock, string(s,"1"), "Q")
            end
            print("call : S", s, " end \n")
        end

        if (mode == "R")
            print("call : R")
            if (k-3 >= 0)
                X1 = FDACA(x, alpha[2:k-1], betta[2:k-2], lenBlock, string(s,"0"), "R")
            end
            if (len-(k+1) >= 0)
                X2 = FDACA(x, alpha[k+1:end], betta[k+1:end], lenBlock, string(s,"1"), "T")
            end
            if (k-4 >= 0)
                X3 = FDACA(x, alpha[2:k-2], betta[2:k-3], lenBlock, string(s,"0"), "Q")
            end
            if (len-(k+2) >= 0)
                X4 = FDACA(x, alpha[k+2:end], betta[k+2:end], lenBlock, string(s,"1"), "R")
            end
            print("call : R", s, " end \n")
        end

        if (mode == "Q")
            print("call : Q")
            if (k-3 >= 0)
                X1 = FDACA(x, alpha[2:k-1], betta[2:k-2], lenBlock, string(s,"0"), "R")
            end
            if ((len-1)-(k+1) >= 0)
                X2 = FDACA(x, alpha[k+1:end-1], betta[k+1:end-1], lenBlock, string(s,"1"), "S")
            end
            if ((k-2)-(2) >= 0)
                X3 = FDACA(x, alpha[2:k-2], betta[2:k-3], lenBlock, string(s,"0"), "Q")
            end
            if ((len-1)-(k+2) >= 0)
                X4 = FDACA(x, alpha[k+2:end-1], betta[k+2:end-1], lenBlock, string(s,"1"), "Q")
            end
            print("call : Q", s, " end \n")
        end

        X0 = X1*X2 + abs(betta[k])*X1*X4 + abs(betta[k])*X2*X3
        return X0
    end
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

#=
function bulgeChasing(tridiagonal)
n = size(tridiagonal)[2]
b = zeros(n,n)
k = n-2
p = 4
N = (n-2)/p
for l = 1:1:N
    s = (l - 1)*p+1
    for j = s:1:(s+p-1)
        q = s + mod(j-1, p) + 1
        kq = q + k - 1
        b[j] =


end=#

A = getGodunovMatrix(16)[1]

#singVals1 = svd(A)[2]

#tridiagonal = lanczosFull(A)
listsAlphAndBetta = lanczos(A)
alpha = Float64.(listsAlphAndBetta[1])
betta = Float64.(listsAlphAndBetta[2])
tridiagonal = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])
#=singVals2 = sqrt.(svd(tridiagonal)[2])

residual = singVals1 - singVals2
#file = open("residualLanczosSingular.txt", "w")
#writedlm(file, residual)
#close(file)
=#
#

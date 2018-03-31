
function getGodunovMatrix(godunov_matrix_dim)

n = godunov_matrix_dim
A = fill(0e0,n,n)  # Матрица A размерности n*n типа BigFloat задается
                       # функцией fill() и заполняется нулевыми элементами.
B = Array(Float64, n) # Вектор B – массив размерности n типа BigFloat.

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

real_result = Array(Float64, n)
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
   	alpha = Array{Float64}(n)

   	### betta - subdiagonal elements
   	betta = Array{Float64}(n-1)

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

#=
function asd(alpha, betta, len)

#len = length(alpha)
#print("\n length(alpha) ->", length(alpha), " \n")
if (len <= 4)
    #show(alpha)
    #a = [alpha[1] betta[1]; betta[1] alpha[2]]
    a = toDense(alpha, betta)
    #show(a)
    show(svd(a)[2])
    #show([alpha[1] betta[1]; betta[2] alpha[2]])
    #push!(listEigvals, eigvals([alpha[1] betta[1]; betta[2] alpha[2]]))
else
    k = div(len, 2)
    #print(k, "  ")

    alpha[k] = alpha[k]-abs(betta[k])
    alpha[k+1] = alpha[k+1]-abs(betta[k])
    #print("start ", k ,"\n")
    asd(alpha[1:k], betta[1:k-1], k)
    #print("Luck! ", k, " end\n")
    asd(alpha[k+1:end], betta[k+1:end], k)
end
end=#

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

A = getGodunovMatrix(1000)[1]

singVals1 = svd(A)[2]

listsAlphAndBetta = lanczos(A)
tridiagonal = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])
singVals2 = sqrt.(svd(A)[2])

residual = singVals1 - singVals2
file = open("residualLanczosSingular.txt", "w")
writedlm(file, residual)

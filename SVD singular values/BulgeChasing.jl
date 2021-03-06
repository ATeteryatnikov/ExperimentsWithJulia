""" Получение матрицы Гивенса """
function givensMatrix(size, a, b, numDiagonalElement=1, blockSize=1)

    G = eye(size)
    r = sqrt(a^2+b^2)
    if (r==0)
        return "Error r=0"
    end
    c = a/r
    s = -b/r
    G[numDiagonalElement, numDiagonalElement] = c
    G[numDiagonalElement, numDiagonalElement+blockSize] = -s
    G[numDiagonalElement+blockSize, numDiagonalElement] = s
    G[numDiagonalElement+blockSize, numDiagonalElement+blockSize] = c

    return G
end

""" Поиск сингулярных чисел трехдиагональной симметричной матрицы A (плотной)"""
function bulgeChasing(A, starterNumElement, repeatCount)

    rowsCount = size(A)[1]

    for iter in 1:1:repeatCount

        GM = givensMatrix(rowsCount, 
                         A[starterNumElement,starterNumElement],
                         A[starterNumElement+1,starterNumElement], starterNumElement, 1)
        A = GM*A*GM'

        for i in starterNumElement:1:(rowsCount-2)
            GM = givensMatrix(rowsCount, A[i+1,i], A[i+2,i], i+1, 1)
            A = GM*A*GM'
        end

        #GM = givensMatrix(rowsCount, A[rowsCount-1,rowsCount-1], A[rowsCount,rowsCount-1], rowsCount-1, 1)
        #A = GM*A*GM'

    end
    A = setZeroElements(A, 1e-10)
    return A
end

function setZeroElements(A, minVal)
    m = size(A)[1]
    for i in 1:1:m-2
        for j in i+2:1:m
            if (abs(A[i,j]) < minVal)
                A[i,j] = 0
            end
            if (abs(A[j,i]) < minVal)
                A[j,i] = 0
            end
        end
    end
    return A
end

function getGVals(a, b)
    r = sqrt(a^2+b^2)
    if (r==0)
        return "Error r==0"
    end
    c = a/r
    s = -b/r
    return [c,s]
end

function sparseBulgeChasing(newAlphas, newBettas, numRepeat, startElement)
  
    N = length(newAlphas)  
for numRep in 1:1:numRepeat
    # значения cos и sin для матрицы Гивенса (1,2).
        (cos0, sin0) = getGVals(newAlphas[1],newBettas[1])
    # после умножения трехдиагональной матрицы слева и справа на транспонированую матрицу Гивенса
    # полученная матрица теряет трехдиагональность из-за одного элемента выходящего за диагонали. 
    # Обозначим выступающий элемент 'P'.
    # Вычислим матричное умножение для симметричной матрицы по формулам полученным в mathcad:
        a1 = newAlphas[1]
        a2 = newAlphas[2]
        b1 = newBettas[1]
        b2 = newBettas[2]

        newAlphas[1] = a1*cos0^2 - 2*b1*cos0*sin0 + a2*sin0^2
        newBettas[1] = b1*cos0^2 - b1*sin0^2 + a1*cos0*sin0 - a2*cos0*sin0
        newAlphas[2] = a2*cos0^2 + 2*b1*cos0*sin0 + a1*sin0^2
        newBettas[2] = b2*cos0
        P = -b2*sin0

        for i in 1:1:length(alphas)-3
        # Для исключения выступающего элемента P вычисляем матрицу Гивенса (i+1,i+2)
            (cosN, sinN) = getGVals(newBettas[i], P)

            b1 = newBettas[i]
            b2 = newBettas[i+1]
            b3 = newBettas[i+2]
            a1 = newAlphas[i]
            a2 = newAlphas[i+1]
            a3 = newAlphas[i+2]
        # Умножаем слева и справа на матрицу Гивенса.
            newBettas[i] = b1*cosN - P*sinN
            newBettas[i+1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
            newBettas[i+2] = b3*cosN
            newAlphas[i+1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
            newAlphas[i+2] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
        # Получаем новый выступающий элемент P на позиции (i+2,i)       
            P = -b3*sinN
        end

    # Последняя итераци исключает P без появления нового

        (cosN, sinN) = getGVals(newBettas[N-2], P)

        b1 = newBettas[N-2]
        b2 = newBettas[N-1]
        a1 = newAlphas[N-2]
        a2 = newAlphas[N-1]
        a3 = newAlphas[N]

        newBettas[N-2] =  b1*cosN - P*sinN
        newBettas[N-1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
        newAlphas[N-1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
        newAlphas[N] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2

    end

    fullMatr = Float64.(toDense(newAlphas, newBettas))

    return [newAlphas, newBettas]
end


##  Демонстрация
n=50

include("Lanczos.jl")
godunovLists = lanczos(getGodunovMatrix(n)[1])
alphas = godunovLists[1]
bettas = godunovLists[2]

Af = Array{Float64}(toDense(copy(alphas), copy(bettas)))

sparseBCRes = sparseBulgeChasing(copy(alphas), copy(bettas),10000,1)

svdValsBC = abs.(sparseBCRes[1])
#originVals = svdvals(Af)
#difference = originVals - svdValsBC


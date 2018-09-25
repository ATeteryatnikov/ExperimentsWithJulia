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

##  Демонстрация метода

alphas = [1 3 0.000007 9999 5 8]
bettas = [2 500 9 1 7]
Af = Array{Float64}(toDense(alphas, bettas))


G1 = givensMatrix(6, Af[1,1], Af[2,1], 1, 1)
A1 = G1*Af

G2 = givensMatrix(6, A1[2,2], A1[3,2], 2, 1)
A2 = G2*G1*Af

G3 = givensMatrix(6, A2[3,3], A2[4,3], 3, 1)
A3 = G3*G2*G1*Af

G4 = givensMatrix(6, A3[4,4], A3[5,4], 4, 1)
A4 = G4*G3*G2*G1*Af

G5 = givensMatrix(6, A4[5,5], A4[6,5], 5, 1)
A5 = G5*G4*G3*G2*G1*Af

G = G5*G4*G3*G2*G1

Agen = copy(Af)
Ggen = eye(6)
Gt=0



for k in 1:1:2000
Ggen = eye(6)
for i in 1:1:(size(Af)[1]-1)
    #print(i)
    Gt = givensMatrix(6, Agen[i,i], Agen[i+1,i], i, 1)
    Agen = Gt*Agen
    Ggen = Gt*Ggen
end

Agen = Agen*Ggen'
end



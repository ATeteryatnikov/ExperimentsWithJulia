try
    include("Lanczos.jl")
    include("DevideAndConqueror.jl")
catch ex
    print("Что-то пошло не так: ", ex)
end

setprecision(200)
function bisection(x, y, alpha, betta, lenBlock)
    k = 0
    x = big(x)
    y = big(y)
    z = big(0.0)
    if (sign(FDACA(x, copy(alpha), copy(betta), lenBlock))==sign(FDACA(y, copy(alpha), copy(betta), lenBlock)))
        return "Знаки совпали"
    end

    while ( (y-x)>1e-15 )
        k = k + 1
        h = (y-x)/2
        z = x + h

        #print("\nX: ", Float64(x),"Z: ", Float64(z),"Y: ", Float64(y))
        fdacaZ = FDACA(z, copy(alpha), copy(betta), lenBlock)
        znakZ = Int(sign(fdacaZ))
        fdacaX = FDACA(x, copy(alpha), copy(betta), lenBlock)
        znakX = Int(sign(fdacaX))
        fdacaY = FDACA(y, copy(alpha), copy(betta), lenBlock)
        znakY = Int(sign(fdacaY))
        if (znakX == znakY)
            print("\nОШИБКА ТУТ\n")
            print("\nX: ", Float64(x)," Z: ", Float64(z)," Y: ", Float64(y))
            print("FDACAx: ", fdacaX,"\n")
            print("FDACAy: ", fdacaY,"\n")
            #print("FDACAz \n", fdacaZ,"\n")
            throw("знаки совпали X=Y (1)")
        end
        #print(znak)
        #print("\n znakX, znakZ, znakY: ", znakX, " " , znakZ, " ", znakY)
        if (znakZ == znakX)
            x = z
        #    print(" X == Z")
        elseif (znakZ == znakY)
            y = z
            #print(" Y == Z")
        end

    end
    #print(k, "\n")
    return z
end

###########################  Check FDACA with Bisection test

#=
x = -1
lenBlock = 4

alphas = [1 3 0.000007 9999 5 8]
bettas = [2 500 9 1 7]
Af = toDense(alphas, bettas)
res1 = det(Af - x*eye(6))
res2 = FDACA(x, copy(alphas), copy(bettas), 4)
diffres = res1 - res2

eigVal = eigvals(Af)[2]
print("\nOne eigen value is: ", eigVal)

d = det(Af - eigVal*eye(6))
print("\n determinant in one eigen value point: ", d)
detFDACA = FDACA(eigVal, copy(alphas), copy(bettas), 4)
print("\n determinant by FDACA in one eigen value point: ", detFDACA)

oneValue = bisection(-1,0.5, alphas, bettas, lenBlock)
print("\n метод бисекций получает: ",oneValue)
fdacaRes = FDACA(oneValue, copy(alphas), bettas, lenBlock)
print("\n Определитель в точке = ", fdacaRes)=#
###########################  Check FDACA


#=
bigFloatMantissa = 100
n = 100
x = sqrt(656.6815092957207)
lenBlock = 4

setprecision(bigFloatMantissa)
A = getGodunovMatrix(n)[1]

listsAlphAndBetta = lanczos(A)
denseTridiagonalMatrix = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])

detTridiagonalMatrix = det(denseTridiagonalMatrix - x*eye(n))
detByFDACA = FDACA(x, listsAlphAndBetta[1], listsAlphAndBetta[2], lenBlock)

diffDet = detTridiagonalMatrix - detByFDACA
=#

x = -1
lenBlock = 4
n = 40

godunov = getGodunovMatrix(n)
alphas = diag(godunov[1])
bettas = diag(godunov[1],1)


Af = toDense(alphas, bettas)
res1 = det(Af - x*eye(n))
res2 = FDACA(x, copy(alphas), copy(bettas), 4)
diffres = res1 - res2

eigVal = eigvals(Af)[2]
print("\nOne eigen value is: ", eigVal)

d = det(Af - eigVal*eye(n))
print("\n determinant in one eigen value point: ", d)
detFDACA = FDACA(eigVal, copy(alphas), copy(bettas), 4)
print("\n determinant by FDACA in one eigen value point: ", detFDACA)

oneValue = bisection(1,1.3, copy(alphas), bettas, lenBlock)
print("\n метод бисекций получает: ",oneValue)
fdacaRes = FDACA(oneValue, copy(alphas), bettas, lenBlock)
print("\n Определитель в точке = ", fdacaRes)

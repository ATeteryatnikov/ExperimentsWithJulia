try
    include("src/Matrices.jl")
    include("Lanczos.jl")
    include("DevideAndConqueror.jl")
    using GenericSVD
catch ex
    print("Что-то пошло не так: ", ex)
end

numTest = 5

# 1)
###########################  Check FDACA Bisection test
if (numTest == 1)
    print("\nТест 1) Ищем сингулярное число матрицы 6х6 методом бисекций и ",
         "\nнаходим значение характеристического многочлена в этой точке")
    setprecision(300)
    lenBlock = 4

    alphas = [1 3 0.000007 9999 5 8]
    bettas = [2 500 9 1 7]
    Af = toDense(alphas, bettas)

    eigVal = svdvals(Af)[2]
    #print("\nOne eigen value is: ", eigVal)

    detInOneEigVal = det(Af - eigVal*eye(6))
    #print("\n determinant in one eigen value point: ", detInOneEigVal)
    detFDACA = FDACA(eigVal, copy(alphas), copy(bettas), 4)
    #print("\n determinant by FDACA in one eigen value point: ", detFDACA)
    print("\n Разность значений  характеристических многочленов",
          "\n полученных разными способами в точке сингулярного числа: ", detFDACA-detInOneEigVal)

    OneEigenValue = bisection(-1,0.5, 1e-20, copy(alphas), bettas, lenBlock)
    print("\n метод бисекций получает сингулярное число: ", OneEigenValue)
    fdacaRes = FDACA(OneEigenValue, copy(alphas), bettas, lenBlock)
    print("\n Характеристический многочлен принимает значение: ", fdacaRes)
end

# 2)
############    Check on Godunov matrix
if (numTest == 2)
    print("\nТест 2) Ищем сингулярное число матрицы Годунова методом бисекций и ",
         "\nнаходим значение характеристического многочлена в этой точке")
    setprecision(500)

    x = -1
    lenBlock = 4
    n = 10

    godunovLists = lanczos(getGodunovMatrix(n)[1])
    alphas = godunovLists[1]
    bettas = godunovLists[2]


    eigVal = svdvals(Af)[2]

    detInOneEigVal = det(Af - eigVal*eye(6))
    #print("\n determinant in one eigen value point: ", detInOneEigVal)
    detFDACA = FDACA(eigVal, copy(alphas), copy(bettas), 4)
    #print("\n determinant by FDACA in one eigen value point: ", detFDACA)
    print("\n Разность значений  характеристических многочленов",
          "\n полученных разными способами в точке сингулярного числа: ", detFDACA-detInOneEigVal)

    OneEigenValue = bisection(-1,0.5, big(1e-20), copy(alphas), bettas, lenBlock)
    print("\n метод бисекций получает сингулярное число: ", OneEigenValue)
    fdacaRes = FDACA(OneEigenValue, copy(alphas), bettas, lenBlock)
    print("\n Характеристический многочлен принимает значение: ", fdacaRes)
end

# 3)
###########################
if (numTest == 3)
    print("\nТест 3) Ищем значение характеристического многочлена в произвольной точке ",
         "\nдля трехдиагональной матрицы полученной из матрицы матрицы Годунова методом Ланцоша")
    setprecision(100)
    x = 20

    bigFloatMantissa = 100
    n = 10
    lenBlock = 4

    setprecision(bigFloatMantissa)
    A = getGodunovMatrix(n)[1]

    listsAlphAndBetta = lanczos(A)
    denseTridiagonalMatrix = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])

    detTridiagonalMatrix = det(denseTridiagonalMatrix - x*eye(n))
    detByFDACA = FDACA(x, listsAlphAndBetta[1], listsAlphAndBetta[2], lenBlock)

    diffDet = detTridiagonalMatrix - detByFDACA
    print("\n Разность значений характеристических многочленов в точке: ", diffDet)
end

# 4)
####### test FDACA signum
if (numTest == 4)
    print("\nТест 4) Ищем разными методами значение характеристического многочлена для матрицы 6х6",
         "\nв сингулярных точках матрицы полученных методом svdvals. Сравниваем знак")
    setprecision(100)

    lenBlock = 4

    alphas = [1 3 0.000007 9999 5 8]
    bettas = [2 500 9 1 7]
    Af = toDense(alphas, bettas)
    singOrigin = svdvals(Af)

    for i in 1:1:6
        detOrigin = det(Af - singOrigin[i]*eye(6))
        #print(Float64(detOrigin), "\n")
        detFDACA = FDACA(singOrigin[i], copy(alphas), copy(bettas), 4)
        #print(Float64(detFDACA), "\n")

        if sign(detOrigin)!=sign(detFDACA)
            print("Проблема на сингулярном числе № ",i,"\n")
        end
    end
end

#########   Check svdValsFinder
if (numTest == 5)

    lenBlock = 4

    alphas = [1 3 0.000007 9999 5 8]
    bettas = [2 500 9 1 7]
    Af = toDense(alphas, bettas)

    sing = svdValsFinder(copy(alphas), copy(bettas), lenBlock, 1e-20)
    singOrigin = svdvals(Af)

end

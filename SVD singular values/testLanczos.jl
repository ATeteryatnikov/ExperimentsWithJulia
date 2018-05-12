try
    include("Lanczos.jl")
catch ex
    print("Что-то пошло не так: ", ex)
end

###########################  Check Lanczos Algorithm
bigFloatMantissa = 500
n = 100

setprecision(bigFloatMantissa)
A = getGodunovMatrix(n)[1]

listsAlphAndBetta = lanczos(A)
denseTridiagonalMatrix = toDense(listsAlphAndBetta[1], listsAlphAndBetta[2])

###     Нахождение сигулярных чисел исходной и треугольной
singValsOriginalMatrix = svdvals(A'*A)
singValsTridiagonalMatrix = svdvals(denseTridiagonalMatrix)

diffSingVals = singValsOriginalMatrix - singValsTridiagonalMatrix
normDiff = norm(diffSingVals)

###     Запись в файлы
path = string("resultLanczosTests/dim_", n, "_matissa_", bigFloatMantissa, ".txt")
fileDifferenceSingVals = open(path, "w")
# Записываем результат в формате Float64
writedlm(fileDifferenceSingVals, Float64.(diffSingVals))
close(fileDifferenceSingVals)

path = string("resultLanczosTests/SingVals_dim_", n, "_matissa_", bigFloatMantissa, ".txt")
SingVals = open(path, "w")
# Записываем результат в формате Float64
writedlm(SingVals, Float64.(sqrt.(singValsTridiagonalMatrix)))
close(SingVals)

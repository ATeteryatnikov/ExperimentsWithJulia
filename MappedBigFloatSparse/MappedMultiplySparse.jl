try
    using JLD
catch
    print("Ошибка при загрузке JLD")
end

""" Имя файла логов (файл создастся сам)"""
logFileName = "/logMultiply.txt"

# Подключаются дополнительные файлы, которые должны лежать в каталоге src,
# который находится в каталоге с исполняемым файлом
path = @__FILE__
path = dirname(path)
include(string(path,"/ReadIndexStructureAndConnectMapped.jl"))


""" Имя файла результата матричного умножения (файл создастся сам)"""
resultFileName = "multiplyResult.jld"

RHS = big.(rand(width))

@showAndLog(string("Умножение матрицы на вектор началось..."))
tic()
result = MulOnVec(file, IndStruct, RHS)
timeMultiply = toc()
@showAndLog(string("Вычислилось за : ", timeMultiply, " сек."))

resfile = jldopen(resultFileName, "w")
resfile["res"] = result
resfile["RHS"] = RHS

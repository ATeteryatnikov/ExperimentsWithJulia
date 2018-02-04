try
    using JLD
catch
    print("Ошибка при загрузке JLD")
end

""" Имя файла логов (файл создастся сам)"""
logFileName = "/logMultiply.txt"

""" Имя файла логов (файл создастся сам)"""
RHSFileName = "rhscopy.txt"

# Подключаются дополнительные файлы, которые должны лежать в каталоге src,
# который находится в каталоге с исполняемым файлом
path = @__FILE__
path = dirname(path)
#include(string(path,"/ReadIndexStructureAndConnectMapped.jl"))
include(string(path,"/src/SparseSolvers.jl"))
include(string(path,"/src/SparseNumberType.jl"))

""" Имя jld файла результата матричного умножения (файл создастся сам)"""
resultFileName = "multiplyResult.jld"

BigFloatFunc = (x)->parse(BigFloat, x)
file = readdlm(RHSFileName, String)
RHS = map(BigFloatFunc, file)
file = 0

xStartPoint = big.(zeros(width))
steps=2

@showAndLog(string("Выполнение МСГ началось..."))
tic()
result = doMCG(fileA, IndStruct,
               fileAtr, JndStruct,
               RHS, xStartPoint,  steps)
timeMultiply = toc()
@showAndLog(string("Выполнилось ", steps, " шагов МСГ за : ", timeMultiply, " сек."))

resfile = jldopen(resultFileName, "w")
resfile["res"] = result
resfile["RHS"] = RHS
close(resfile)

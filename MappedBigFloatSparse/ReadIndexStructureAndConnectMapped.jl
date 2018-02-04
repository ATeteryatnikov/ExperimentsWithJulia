try
    using JLD
catch
    print("Ошибка при загрузке JLD")
end


""" Имя файла логов (файл создастся сам)"""
logFileName = "/logRead.txt"

# Подключаются дополнительные файлы, которые должны лежать в каталоге src,
# который находится в каталоге с исполняемым файлом
path = @__FILE__
path = dirname(path)
include(string(path,"/src/LogGenerator.jl"))
include(string(path,"/src/SparseNumberType.jl"))


""" Имя файла со значениями матрицы  (файл создастся сам)"""
nameJlDFileWithValues = "PartitionalJLD.jld"

""" Имя файла с индексной структурой  (файл создастся сам)"""
nameJLDIndexStructure = "indexStructure.jld"

""" Имя файла с индексной структурой для транспонированной матрицы (файл создастся сам)"""
nameJLDIndexStructureTr = "indexStructureTr.jld"

""" Имя файла со значениями матрицы для транспонированной матрицы (файл создастся сам)"""
nameJlDFileWithValuesTr = "PartitionalJLDTr.jld"

@showAndLog("\n Start!-------------------------")

@showAndLog(string("Подключается JLD файл с данными BigFloat с параметром mmaparrays..."))
tic()
fileA = jldopen(nameJlDFileWithValues, "r", mmaparrays=true)
timeConnectionBigFloat = toc()
@showAndLog(string("JLD файл подключился за ", timeConnectionBigFloat, " сек."))

@showAndLog(string("Считывается JLD файл структуры..."))
tic()
indexFile = jldopen(nameJLDIndexStructure, "r")
IndStruct = read( indexFile , "Ind")
height = read( indexFile , "height")
width = read( indexFile , "width")
close(indexFile)
timeReadStructure = toc()
@showAndLog(string("Структура считана за ", timeReadStructure, " сек."))

####################### Для Транспонированной матрицы

@showAndLog(string("Подключается JLD файл с данными BigFloat для траспонированной матрицы
 с параметром mmaparrays..."))
tic()
fileAtr = jldopen(nameJlDFileWithValuesTr, "r", mmaparrays=true)
timeConnectionBigFloatTr = toc()
@showAndLog(string("JLD файл для траспонированной матрицы подключился за ", timeConnectionBigFloatTr, " сек."))

@showAndLog(string("Считывается JLD файл структуры для траспонированной матрицы..."))
tic()
indexFile = jldopen(nameJLDIndexStructureTr, "r")
JndStruct = read( indexFile , "Jnd")
close(indexFile)
timeReadStructureTr = toc()
@showAndLog(string("Структура для траспонированной матрицы считана за ", timeReadStructureTr, " сек."))

@showAndLog(string("Подключены jld файлы: ",nameJlDFileWithValues , ", ",nameJlDFileWithValuesTr))
@showAndLog(string("На чтение в переменные: fileA, fileAtr "))
@showAndLog(string("Считаны индексные структуры в переменные: IndStruct, JndStruct" ))

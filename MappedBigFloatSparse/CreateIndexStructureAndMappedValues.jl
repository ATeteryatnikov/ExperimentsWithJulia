try
    using JLD
catch
    print("Ошибка при загрузке JLD")
end

""" Имя файла с исходной матрицей (файл должен уже существовать!) """
fileName = "/matrix.txt"

""" Имя файла логов (файл создастся сам)"""
logFileName = "/logCreate.txt"

""" Имя файла со значениями матрицы  (файл создастся сам)"""
nameJlDFileWithValues = "PartitionalJLD.jld"

""" Имя файла со значениями матрицы для транспонированной матрицы (файл создастся сам)"""
nameJlDFileWithValuesTr = "PartitionalJLDTr.jld"

""" Имя файла с индексной структурой  (файл создастся сам)"""
nameJLDIndexStructure = "indexStructure.jld"

""" Имя файла с индексной структурой для транспонированной матрицы (файл создастся сам)"""
nameJLDIndexStructureTr = "indexStructureTr.jld"

# Подключаются дополнительные файлы, которые должны лежать в каталоге src,
# который находится в каталоге с исполняемым файлом
path = @__FILE__
path = dirname(path)
include(string(path,"/src/LogGenerator.jl"))
include(string(path,"/src/SparseNumberType.jl"))


tic()
@showAndLog(string("Чтение файла: ", fileName[2:end], "..."))
TextSourceFile = readdlm(string(path, fileName), String)
timeReadSparsefile = toc()
@showAndLog(string("Считан за: ", timeReadSparsefile, " сек."))

@showAndLog(string("Подготовка данных..."))
tic()
# Функция преобразования текста в число
IntFunc = (x)->parse(Int64, x)+1
# Ind -Все номера строк из фалаю. Jnd - Все номера столбцов в файлк.
Ind = map(IntFunc, TextSourceFile[:,1])
Jnd = map(IntFunc, TextSourceFile[:,2])

# Удаляем из памяти первые 2 столбца TextSourceFile
TextSourceFile = TextSourceFile[:,3]

# Количество строк и столбцов в матрице
height = maximum(Ind)
width = maximum(Jnd)

timeConvert = toc()
@showAndLog(string("Подготовка заняла: ", timeConvert, " сек."))

@showAndLog(string("Формируется структура индексов для матрицы..."))
tic()
IndStruct = createIndexStructure(Ind, Jnd)
timeInd = toc()
@showAndLog(string("Структура индексов сформирована за ", timeInd, " сек."))

@showAndLog(string("Записываем структуру индексов в файл..."))
tic()
file = jldopen(nameJLDIndexStructure, "w")
file["Ind"] = IndStruct
file["width"] = width
file["height"] = height
close(file)
timeWriteInd = toc()
@showAndLog(string("Структура для индексов записана в файл за ", timeWriteInd, " сек."))

@showAndLog(string("Записываем столбец значений типа BigFloat в JLD..."))
tic()
writeInOneJLDGroup(nameJlDFileWithValues, TextSourceFile, IndStruct)
timeWriteBigFloat = toc()
@showAndLog(string("Значения записаны в файл за ", timeWriteBigFloat, " сек."))

################ Для транспонированной матрицы #####################

@showAndLog(string("Формируется структура индексов для транспонированной матрицы..."))
tic()
JndStruct = createIndexStructure(Jnd, Ind)
timeJnd = toc()
@showAndLog(string("Структура индексов сформирована за ", timeJnd, " сек."))

@showAndLog(string("Записываем структуру индексов для транспонированной матрицы в файл..."))
tic()
file = jldopen(nameJLDIndexStructureTr, "w")
file["Jnd"] = JndStruct
file["width"] = width
file["height"] = height
close(file)
timeWriteInd = toc()
@showAndLog(string("Структура индексов для транспонированной матрицы записана в файл за ", timeWriteInd, " сек."))

@showAndLog(string("Записываем столбец значений типа BigFloat для разряженной матрицы в JLD..."))
tic()
writeInOneJLDGroup(nameJlDFileWithValuesTr, TextSourceFile, JndStruct)
timeWriteBigFloat = toc()
@showAndLog(string("Значения записаны в файл за ", timeWriteBigFloat, " сек."))


@showAndLog(string("Созданы jld файлы с значениями BigFloat: ",nameJlDFileWithValues , ", ",nameJlDFileWithValuesTr))
@showAndLog(string("Созданы jld файлы с структурами: ", nameJLDIndexStructure, ", ", nameJLDIndexStructureTr ))


"""макрос выводит информацию на дисплей и в лог-файл. предполагает, что logFileName и path определены до использования макроса"""
macro showAndLog(str)
    logfile = open(string(path, logFileName), "a+")
    str = string(eval(str), "\n")
    write(logfile, str)
    showall(str)
    close(logfile)
end

try
    type myType{}
        i::Int
        jArray::Array{Int}
        rowNumberArray::Array{Int}
    end
catch
    printf("myType уже подключен!")
end

#каждому уникальному значению из Ind сопостовляет массив номеров вхождений этих значений в массив Jnd
function writeIndexStructure(Ind::Array{Int}, Jnd::Array{Int})

    uniqueInd = sort!(unique(Ind))
    n = length(uniqueInd)
    ar = Array{myType}(n)
    for i in 1:1:n
        num = uniqueInd[i]
        rowNumberArray = findin(Ind, uniqueInd[i])
        jArray = Jnd[rowNumberArray]
        ar[i] = myType(num, rowNumberArray, jArray)
    end
    return ar
end

function multipleOnDenseVector(Ind, rhs)

    rowsA = length(Ind)
    result = big(zeros(rowsA))
    for i in 1:1:rowsA
    	#Перебираем все ненулевые столбцы
    	#Ind[i].jArray - массив ненулевых столбцов, что Array[i,j]!=0
        #print(i, " ")
    	   result[i] = (rhs[Ind[i].jArray]'*map((x)->BigFloatArray["val"][x], Ind[i].rowNumberArray))[1]
    end
end

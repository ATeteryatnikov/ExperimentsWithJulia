""" функция возвращает массив элементы которого:
1) матрица Годунова (левая часть СЛАУ Годунова)
2) вектор - правая часть СЛАУ Годунова
3) вектор - решение СЛАУ Годунова
Источник: Годунов С.К. Решение систем линейных уравнений - Новосибирск: Наука, 1980"""
function getGodunovMatrix(n)

    # создается нулевая Матрица A размерности n*n типа BigFloat
    A = fill(big(0e0), n, n)
    # Вектор B – массив размерности n типа BigFloat.
    B = Array{BigFloat}(n)

    # При заполнении элементов массива, необходимо сначала
    # преобразовать целое число в тип BigFloat, и только
    # после этого, делить его на другое число.
    diagonalElement = big(7)/5
    upperDiagonalElement = big(11)/3

    # заполнение ненулевых элементтов матрицы A
    for i in 1:1:(n-1)
        A[i,i] = diagonalElement
        A[i,i+1] = upperDiagonalElement
    end
    A[n,n] = diagonalElement

    # Заполнение вектора B по формулам
    for i in 1:1:(n-2)
        B[i] = big(152*i+118)
        B[i] = B[i]/(15*(2*i+1)*(2*i+3))
    end
    B[n-1] = big(152*n-34)/(15*(2*n-1)*(2*n+1))
    B[n] = big(7)/(5*(2*n+1))

    # Вычисления решения системы Ax=B
    # по формулам для нахождения решения СЛАУ Годенова
    real_result = Array{BigFloat}(n)

    for i in 1:1:(n-2)
        real_result[i] = 1/(2*i+1)
    end
    real_result[n-1] = 1/(2*n-1)
    real_result[n] = 1/(2*n+1)

    return [A, B, real_result]
end

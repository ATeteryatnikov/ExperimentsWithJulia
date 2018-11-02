""" функция возвращает косинус и синус - элементы матрицы Гивенса  """
function getGVals(a, b)
    r = sqrt(a^2+b^2)
    if (r==0)
        return "Error r == 0"
    end
    c = a/r
    s = -b/r
    return [c,s]
end

function sparseBulgeChasing(newAlphas, newBettas, numRepeat, startElement)
  
    N = length(newAlphas)  
    for numRep in 1:1:numRepeat
    # значения cos и sin для матрицы Гивенса (1,2).
        (cos0, sin0) = getGVals(newAlphas[1],newBettas[1])
    # после умножения трехдиагональной матрицы слева и справа на транспонированую матрицу Гивенса
    # полученная матрица теряет трехдиагональность из-за одного элемента выходящего за диагонали. 
    # Обозначим выступающий элемент 'P'.
    # Вычислим матричное умножение для симметричной матрицы по формулам полученным в mathcad:
        a1 = newAlphas[1]
        a2 = newAlphas[2]
        b1 = newBettas[1]
        b2 = newBettas[2]

        newAlphas[1] = a1*cos0^2 - 2*b1*cos0*sin0 + a2*sin0^2
        newBettas[1] = b1*cos0^2 - b1*sin0^2 + a1*cos0*sin0 - a2*cos0*sin0
        newAlphas[2] = a2*cos0^2 + 2*b1*cos0*sin0 + a1*sin0^2
        newBettas[2] = b2*cos0
        P = -b2*sin0

        for i in 1:1:length(alphas)-3
        # Для исключения выступающего элемента P вычисляем матрицу Гивенса (i+1,i+2)
            (cosN, sinN) = getGVals(newBettas[i], P)

            b1 = newBettas[i]
            b2 = newBettas[i+1]
            b3 = newBettas[i+2]
            a1 = newAlphas[i]
            a2 = newAlphas[i+1]
            a3 = newAlphas[i+2]
        # Умножаем слева и справа на матрицу Гивенса.
            newBettas[i] = b1*cosN - P*sinN
            newBettas[i+1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
            newBettas[i+2] = b3*cosN
            newAlphas[i+1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
            newAlphas[i+2] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
        # Получаем новый выступающий элемент P на позиции (i+2,i)       
            P = -b3*sinN
        end

    # Последняя итераци исключает P без появления нового

        (cosN, sinN) = getGVals(newBettas[N-2], P)

        b1 = newBettas[N-2]
        b2 = newBettas[N-1]
        a1 = newAlphas[N-2]
        a2 = newAlphas[N-1]
        a3 = newAlphas[N]

        newBettas[N-2] =  b1*cosN - P*sinN
        newBettas[N-1] = b2*cosN^2 - b2*sinN^2 + a2*cosN*sinN - a3*cosN*sinN
        newAlphas[N-1] = a2*cosN^2 - 2*b2*cosN*sinN + a3*sinN^2
        newAlphas[N] = a3*cosN^2 + 2*b2*cosN*sinN + a2*sinN^2
    end

    return [newAlphas, newBettas]
end


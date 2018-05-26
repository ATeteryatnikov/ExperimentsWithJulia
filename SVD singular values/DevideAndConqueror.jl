""" Вычисление полинома в точке x по схеме  4.6 из статьи V. Rokhlin
A fast divide-and-conquer algorithm for computing the spectra of real symmetric tridiagonal matrices """
function FDACA(x, alpha, betta, lenBlock, s="", mode="T")

    level = length(s)
    len = length(alpha)
    #print("Level: ", s, " Mode : ", mode, " lenAlpha: ", len, " lenBetta: ", length(betta), "\n")
    #print("WRONG:  ", (len<=lenBlock))
    if (len <= lenBlock)
        #print("LENGTH RIGHT: ", len,"!\n")
        matr = toDense(alpha, betta)

        return det(matr - x*eye(len))
    else
        k = div(len, 2)
        #print("\n",k, "is K \n  ")

        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k+1] = alpha[k+1] - abs(betta[k])

        X1 = X2 = X3 = X4 = 1

        if (mode == "T")
            #print("call : T \n")
            if (k-1 >= 0)
                X1 = FDACA(x, copy(alpha[1:k]), copy(betta[1:k-1]), lenBlock, string(s,"0"), "T")
                #matr = toDense(alpha[1:k], betta[1:k-1])
                #X1 = det(matr - x*eye(k))
                #print("X1: ", X1, "\n")
                #show(matr)
                #print("\n")
            end
            if (len-(k+1) >= 0)
                X2 = FDACA(x, copy(alpha[k+1:end]), copy(betta[k+1:end]), lenBlock, string(s,"1"), "T")
                #matr = toDense(alpha[k+1:end], betta[k+1:end])
                #X2 = det(matr - x*eye(len-k))

                #print("X2: ", X2, "\n")
                #return matr
                #show(matr)
                #print("\n")
                #print("len alpha = ", length(alpha[k+1:end]), "\n")
            end
            if (k-2 >= 0)
                X3 = FDACA(x, copy(alpha[1:k-1]), copy(betta[1:k-2]), lenBlock, string(s,"0"), "S")
                #matr = toDense(alpha[1:k-1], betta[1:k-2])
                #X3 = det(matr - x*eye(k-1))

                #print("X3: ", X3, "\n")
                #print("len alpha = ", length(alpha[1:k-1]), "\n")
                #show(matr)
                #print("\n")
            end
            if (len-(k+2) >= 0)
                X4 = FDACA(x, copy(alpha[k+2:end]), copy(betta[k+2:end]), lenBlock, string(s,"1"), "R")
                #matr = toDense(alpha[k+2:len], betta[k+2:len-1])
                #X4 = det(matr - x*eye(len-k-1))

                #print("X4: ", X4, "\n")
                #print("len alpha = ", length(alpha[k+2:end]), "\n")
                #show(matr)
                #print("\n")
            end
        end

        if (mode == "S")
            #print("call : S")
            if (k-1 >= 0)
                X1 = FDACA(x, copy(alpha[1:k]), copy(betta[1:k-1]), lenBlock, string(s,"0"), "T")
            end
            if ((len-1)-(k+1) >= 0)
                X2 = FDACA(x, copy(alpha[k+1:end-1]), copy(betta[k+1:end-1]), lenBlock, string(s,"1"), "S")
            end
            if (k-2 >= 0)
                X3 = FDACA(x, copy(alpha[1:k-1]), copy(betta[1:k-2]), lenBlock, string(s,"0"), "S")
            end
            if ((len-1)-(k+2) >= 0)
                X4 = FDACA(x, copy(alpha[k+2:end-1]), copy(betta[k+2:end-1]), lenBlock, string(s,"1"), "Q")
            end
            #print("call : S", s, " end \n")
        end

        if (mode == "R")
            #print("call : R")
            if (k-3 >= 0)
                X1 = FDACA(x, copy(alpha[2:k-1]), copy(betta[2:k-2]), lenBlock, string(s,"0"), "R")
            end
            if (len-(k+1) >= 0)
                X2 = FDACA(x, copy(alpha[k+1:end]), copy(betta[k+1:end]), lenBlock, string(s,"1"), "T")
            end
            if (k-4 >= 0)
                X3 = FDACA(x, copy(alpha[2:k-2]), copy(betta[2:k-3]), lenBlock, string(s,"0"), "Q")
            end
            if (len-(k+2) >= 0)
                X4 = FDACA(x, copy(alpha[k+2:end]), copy(betta[k+2:end]), lenBlock, string(s,"1"), "R")
            end
            #print("call : R", s, " end \n")
        end

        if (mode == "Q")
            #print("call : Q")
            if (k-3 >= 0)
                X1 = FDACA(x, copy(alpha[2:k-1]), copy(betta[2:k-2]), lenBlock, string(s,"0"), "R")
            end
            if ((len-1)-(k+1) >= 0)
                X2 = FDACA(x, copy(alpha[k+1:end-1]), copy(betta[k+1:end-1]), lenBlock, string(s,"1"), "S")
            end
            if ((k-2)-(2) >= 0)
                X3 = FDACA(x, copy(alpha[2:k-2]), copy(betta[2:k-3]), lenBlock, string(s,"0"), "Q")
            end
            if ((len-1)-(k+2) >= 0)
                X4 = FDACA(x, copy(alpha[k+2:end-1]), copy(betta[k+2:end-1]), lenBlock, string(s,"1"), "Q")
            end
            #print("call : Q", s, " end \n")
        end
        #print("lev: ", mode, s, " x1", X1, " x2", X2, " x3", X3, " x4", X4, "\n")
        X0 = X1*X2 + abs(betta[k])*X1*X4 + abs(betta[k])*X2*X3
        #print("\n X0 is\n x1: ", X1, " *  x2 ", X2, " x3 ", X3, " x4 ", X4, " bettak ", betta[k], "\n" )
        #print("call : ",mode, s, " end \n")

        #=sig = sign(X2)*sign(X3)*sign(X1)*sign(X4)
        #print(sig, "\n")
        psi = -log(abs(betta[k])) - log(abs(exp(X4-X2) + sig*exp(X3-X1)))
        if (sign(psi)==sign(X0))
            print("Совпали\n")
        end=#
        return X0
    end
end

""" Поиск сингулярных чисел между двумя заданными числами """
function bisection(x, y, prec, alpha, betta, lenBlock)
    k = 0

    x = big(x)
    y = big(y)
    z = big(0.0)

    if (sign(FDACA(x, copy(alpha), copy(betta), lenBlock))==sign(FDACA(y, copy(alpha), copy(betta), lenBlock)))
        return 0
    end

    while ( (y-x)>prec )
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

""" поиск сингулярных чисел между сингулярными числами дочерных матриц (Т0 и Т1) """
function svdValsFinder(alpha, betta, lenBlock, prec)

    len = length(alpha)

    if (len <= lenBlock)

        matr = toDense(alpha, betta)
        return svdvals(matr)
    else
        k = div(len, 2)

        alpha[k] = alpha[k] - abs(betta[k])
        alpha[k+1] = alpha[k+1] - abs(betta[k])

        if (k-1 >= 0)
            singValsT0 = svdValsFinder(copy(alpha[1:k]), copy(betta[1:k-1]), lenBlock, prec)
        end
        if (len-(k+1) >= 0)
            singValsT1 = svdValsFinder(copy(alpha[k+1:end]), copy(betta[k+1:end]), lenBlock, prec)

        end

        unionSingVals = sort(union(singValsT0, singValsT1), rev=true)
        #return unionSingVals
        singVals = big.(zeros(len))
        #return unionSingVals
        for i in 1:1:len-1
            singVals[i] = bisection(unionSingVals[i+1], unionSingVals[i], prec, copy(alpha), copy(betta), lenBlock)
            print(string(i,": ",Float64(singVals[i]),"\n"))
        end
        singVals[len] = bisection(unionSingVals[len], unionSingVals[len]+2*betta[k], prec, copy(alpha), copy(betta), lenBlock)
        print(string(len,": ",Float64(singVals[len]),"\n"))
        return singVals

    end


end

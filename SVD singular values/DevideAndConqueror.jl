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
        return X0
    end
end

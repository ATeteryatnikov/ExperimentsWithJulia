
include("c:/projects/experiments with Julia/ExperimentsWithJulia/src/Logs.jl")

nameexp = "singValsTest" # Имя эксперимента - добавка спереди к выходным файлам

# setprecision(1000)

matrix = readdlm("c://projects//Matricies//matrix.txt", String)

I = map( (x) -> parse(Int, x), matrix[:,1]) + 1
println("\n\n I: ", typeof(I))
println(size(I))

J = map( (x) -> parse(Int, x), matrix[:,2]) + 1
println("\n\n J: ", typeof(J))
println(size(J))

V = map( (x) -> parse(BigFloat, x), matrix[:,3])
println("\n\n V: ", typeof(V))
println(size(V))


A=sparse(I, J, V)
println("\n\n A: ", typeof(A))
println(size(A))

I = J = V = 0

#=
S = svdvals(A)

println("\n\n\n  sing val = \n", S)
#println("sing val = ", convert(Array{Float64,1},S))

quit()

b = readdlm("rhscopy.txt", Float64)[:,1];
println("\n\n b: ", typeof(b))
println(size(b))
# println(b)

x = readdlm("solution.it1.txt", Float64)[:,1];
println("\n\n x: ", typeof(x))
println(size(x))
# println(x)

nev0=b-A*x

# println("b = ", convert(Array{Float64,1},b))

tic()  #запуск таймера


x1=A\b
nev1=b-A*x1

timer = toc()			# Завершение работы таймера.
println("\n\nsolution timein sec:  ", timer)          # Запись времени работы программы в
                                                  # секундах
flush(STDOUT)

dx=maxabs(x1-x)

#println("\n\n\n dx = ", dx)

#println("\n\n\n x1 = ", x1)

# writedlm("$nameexp.nevyazka0.gnuplot1", nev0 #.*1e6, "# nevyazka (mm)"
#          );

# writedlm("$nameexp.nevyazka1.gnuplot1", nev1 #.*1e6, "# nevyazka (mm)"
#          );

println("\n\nsumabs2(nev0) = ", sumabs2(nev0))
println("\n\nsumabs2(nev1) = ", sumabs2(nev1))

#Pkg.add("GenericSVD")
using GenericSVD

tic()  #запуск таймера

bA = svdfact(A)

timer = toc()			# Завершение работы таймера.
println("time elapsed in sec:  ", timer)          # Запись времени работы программы в
                                                  # секундах

U=bA[:U]
Vt=bA[:Vt]
S=bA[:S]

println("\n\n\n  sing val = \n", S)
#println("sing val = ", convert(Array{Float64,1},S))

A1 = U*diagm(S)*Vt

dA=maximum(abs.(A1-A))
println("dA = ",dA)

x2 = Vt'*diagm(big(1)./S)*U' * b

dx12=maximum(x1-x2)
dx2=maximum(x2-x)

writedlm("$nameexp.x1-x2-dx12.gnuplot1", [x1 x2 dx12] #.*1e6, "# nevyazka (mm)"
         );

writedlm("$nameexp.x-x2-dx2.gnuplot1", [x x2 dx2] #.*1e6, "# nevyazka (mm)"
         );



quit()
=#

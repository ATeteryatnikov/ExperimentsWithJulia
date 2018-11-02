include("Lanczos.jl")
include("SparseBulgeChasing.jl")

for n in 500:500:500
	for mantissa in 100:100:100
		for numRepeated in 100:100:100

			setprecision(mantissa)
			godunovLists = lanczos(getGodunovMatrix(n)[1])
			alphas = godunovLists[1]
			bettas = godunovLists[2]
			fullMatr = toDense(alphas, bettas)
			
			tic()
			# SparseBulgeChasingResult
			SBCRes = sparseBulgeChasing(copy(alphas), copy(bettas), numRepeated, 1)
			timerSBC = toc()

			tic()
			originVals = svdvals(fullMatr)
			denseTimer = toc()

			difference = originVals - abs.(SBCRes[1])
			normDifference = norm(difference)

			# Запись в файлы
			pathDiffSingVals = string("resultBulgeChasingTests/diffSingVals__dim_", n, "_matissa_", mantissa, "_rep_", numRepeated, ".txt")
			fileDifferenceSingVals = open(pathDiffSingVals, "w")
			writedlm(fileDifferenceSingVals, split(string(difference)[10:end-1],","))
			close(fileDifferenceSingVals)

			pathNormDiffSingVals = string("resultBulgeChasingTests/normDiffSingVals__dim_", n, "_matissa_", mantissa, "_rep_", numRepeated, ".txt")
			fileNormDifferenceSingVals = open(pathNormDiffSingVals, "w")
			write(fileNormDifferenceSingVals, string(normDifference))
			close(fileNormDifferenceSingVals)

			pathSBCTime = string("resultBulgeChasingTests/SBCTime__dim_", n, "_matissa_", mantissa, "_rep_", numRepeated, ".txt")
			fileSBCTime = open(pathSBCTime, "w")
			writedlm(fileSBCTime, timerSBC)
			close(fileSBCTime)

			pathDenseTime = string("resultBulgeChasingTests/DenseTime__dim_", n, "_matissa_", mantissa, "_rep_", numRepeated, ".txt")
			fileDenseTime = open(pathDenseTime, "w")
			writedlm(fileDenseTime, denseTimer)
			close(fileDenseTime)
		end
	end
end
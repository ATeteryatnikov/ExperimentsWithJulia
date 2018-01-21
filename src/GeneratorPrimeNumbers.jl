""" (firstPrimeNumber, IterationsLimit = 1, UseTimer = false)\n
Generate big prime number per iterations starts with 'firstPrimeNumber' \n
Default - do 1 iteration\n
Algorith from Osipov N.N."""
function GenerateBigPrimeNumber(firstPrimeNumber, IterationsLimit = 1, UseTimer = false)

	N = firstPrimeNumber
	a = 2
	q = convert(BigInt, N)
	k = 1
	tic()
	while(k <= IterationsLimit)


		alpha = 1:div(q,2)
		R = BigInt(2)*rand(alpha)
		N = q^2*R + 1

		while !(powermod(a, R, N) != 1 && powermod(a, N-1, N) == 1)

			R = 2 * rand(alpha)
			N = q^2 * R + 1
		end

		k = k + 1
		q = N

	end
	UseTimer ? toc():0
	return N
end

cd(@__DIR__)
#using Pkg
#Pkg.activate(".")

using Distributions
using StaticArrays
using Random
using DelimitedFiles
using SLEEF
using JLD2
using FFTW
using Printf


#Random.seed!(parse(Int,ARGS[1]))

### Lattice Size
const L = 8 #parse(Int,ARGS[2])
### end

ξ = Normal(0.0f0, 1.0f0)

### Parameters
const λ = 4.0f0
const Γ = 1.0f0
const T = 1.0f0
### end

### Numerics
const Δt = 0.04f0/Γ
### end

const Rate = Float32(sqrt(2.0*Δt*Γ))

### Nearest neighbors
NNp_a = zeros(Int16,L)
NNm_a = zeros(Int16,L)
for i in 1:L
	NNp_a[i] = i + 1
	NNm_a[i] = i - 1
end
NNp_a[L] = 1
NNm_a[1] = L

# for performance
NNp=@SVector [NNp_a[i] for i in 1:L]
NNm=@SVector [NNm_a[i] for i in 1:L]
###

function hotstart(n)
	rand(ξ, (n, n, n))
end

function ΔH(m², ϕ, ϕt, x)
	@inbounds ϕold = ϕ[x[1],x[2],x[3]]
	Δϕ = ϕt - ϕold
	Δϕ2= ϕt^2 - ϕold^2

	## nearest neighbours
	@inbounds Σnn = ϕ[ NNp[x[1]] , x[2], x[3] ] + ϕ[x[1], NNp[x[2]] ,x[3] ] + ϕ[x[1],x[2], NNp[x[3]] ]
	@inbounds Σnn = Σnn + ϕ[ NNm[x[1]] ,x[2],x[3]] + ϕ[x[1],  NNm[x[2]] ,x[3]] + ϕ[x[1],x[2],  NNm[x[3]] ]

	kinet_term =  3.0f0 * Δϕ2 - Δϕ * Σnn
	poten_term =  0.5f0 * m² * Δϕ2 + 0.25f0 * λ * (ϕt^4-ϕold^4)

	kinet_term + poten_term
end


function MCstep(m², ϕ, x)
	r = rand(Float32)
	@inbounds ϕt = ϕ[x[1],x[2],x[3]] + Rate*rand(ξ)
	p = min(1.0f0,SLEEF.exp(-ΔH(m², ϕ, ϕt, x)))
	if (r<=p)
		@inbounds ϕ[x[1],x[2],x[3]] = ϕt
	end
end


function sweep(m², ϕ, L)
	Threads.@threads for i in 1:L
		for j in 1:L
			for k in 1:L
				if (i+j+k)%2 == 0
					MCstep(m², ϕ, (i,j,k))
				end
			end
		end
	end

	Threads.@threads for i in 1:L
		for j in 1:L
			for k in 1:L
				if (i+j+k)%2 !=0
					MCstep(m², ϕ, (i,j,k))
				end
			end
		end
	end
end

function op(ϕ, L)
	ϕk = fft(ϕ)
	average = ϕk[1,1,1]/L^3
	(real(average),ϕk[:,1,1])
end


function thermalize(m², ϕ, L,  N=10000)
	for t in 1:N
		sweep(m², ϕ, L)
	end
end


###
###
### ************************************** Main ************************************
###
###




## Thermalization study
ϕ=hotstart(L)

# near critical value
m² = -2.285

thermalize(m²,ϕ, L, 100*L^2)

maxt = 10000*L^2

open("output_$L.dat","w") do io 
	for i in 0:maxt
		(M,ϕk) = op(ϕ, L)
		Printf.@printf(io, "%i %f", i, M)
		for kx in 1:L
			Printf.@printf(io, " %f %f", real(ϕk[kx]), imag(ϕk[kx]))
		end 
		Printf.@printf(io, "\n")
		sweep(m², ϕ, L)
	end
end

##

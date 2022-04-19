cd(@__DIR__)

using Plots
using DelimitedFiles
using LaTeXStrings


function autocor_loc_2(x, beg, max, n=2)
	C = zeros(Float64,max+1)
	N = zeros(Int64,max+1)
	Threads.@threads for tau in 0:max
		for i in beg:length(x)-max
			j = i + tau
			@inbounds @fastmath  C[tau+1] = C[tau+1] +  (x[i]*x[j])^n
			@inbounds @fastmath  N[tau+1] = N[tau+1] + 1
		end
	end
	(collect(1:max+1),  C ./ N)
end


df_16=readdlm("output_16.dat",' ')
df_8=readdlm("output_8.dat",' ')

(t_16,c_16) = autocor_loc_2(df_16[:,2], 1, 16^2*100, 1)
(t_8,c_8) = autocor_loc_2(df_8[:,2], 1, 8^2*100, 1)


plot(t_8/8^2,c_8/c_8[1],label=L"L=8",xlabel = L"t/L^2")
plot!(t_16/16^2,c_16/c_16[1],label=L"L=16",ylabel = L"C(t)")

savefig("c.pdf")

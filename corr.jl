cd(@__DIR__)

using Plots
using DelimitedFiles
using LaTeXStrings


function autocor_loc_2(x, beg, max, n=2)
	C = zeros(Complex{Float64},max+1)
	N = zeros(Int64,max+1)
	Threads.@threads for tau in 0:max
		for i in beg:length(x)-max
			j = i + tau
			@inbounds @fastmath  C[tau+1] = C[tau+1] +  (x[i]*conj(x[j]))^n
			@inbounds @fastmath  N[tau+1] = N[tau+1] + 1
		end
	end
	(collect(1:max+1),  C ./ N)
end


df_16=readdlm("output_16.dat",' ')
#df_8=readdlm("output_8.dat",' ')

c_16=Any[]

(t_16,tmp) = autocor_loc_2(df_16[:,2], 1, 16^2, 1)

t_16.=t_16.-1.0

for i in 0:5
	local tmp 
	(tmp_t,tmp) = autocor_loc_2(df_16[:,3+2*i].+df_16[:,4+2*i].*1.0im, 1, 16^2, 1)
	push!(c_16,real.(tmp))
end

#(t_8,c_8) = autocor_loc_2(df_8[:,2], 1, 8^2*100, 1)


#plot(t_8/8^2,c_8/c_8[1],label=L"L=8",xlabel = L"t/L^2")
scatter(20*t_16/16^2,c_16[1]/c_16[1][1],label=L"L=16, C(t,k=0)")
scatter!(20*t_16/16^2,c_16[2]/c_16[2][1],label=L"L=16, C(t,k=2\pi/L)")
scatter!((2.0^2)*20*t_16/16^2,c_16[3]/c_16[3][1],label=L"L=16, C(t,k=4\pi/L)")
scatter!(3^2*20*t_16/16^2,c_16[4]/c_16[4][1],label=L"L=16, C(t,k=6\pi/L)")
#plot!(100*t_16/16^2,c_16[5]/c_16[5][1],label=L"L=16",ylabel = L"C(t)")
scatter!(xlim=(0,10),ylabel = L"C(t)",xlabel=L"t k^z/L^z")



savefig("c.pdf")

### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 017e5156-b4f2-4c90-9bf4-352c7d877ed1
begin
	using Pkg
	Pkg.add(path="/Users/ec627/.julia/dev/Kuramodel")
	using Kuramodel
	greet()
	println("This loads the local Kuramodel Package and keeps it updated.")
end

# ╔═╡ 15f1cffb-e9b8-4abb-892e-48b4e7bea33c
begin
	println("This just calls the packages I am using")
	using Symbolics
	using PlutoUI
	using Distributions
	using Printf
	using Plots
	using Random
	using LinearAlgebra
	using BenchmarkTools
	using Graphs
	using GraphPlot
	using FileIO
	using JLD2
	using LaTeXStrings
	using StatsBase
	using EmpiricalDistributions
end

# ╔═╡ e5277ce7-db25-4a0e-b2fe-f150a0a5ca3b
PlutoUI.TableOfContents()

# ╔═╡ 0f884c5e-a00e-4dda-95ce-89d3d85cf94c
md"""
# Tutorial 1

### Experiment with small network

A good idea would be to start with a system of two oscillators. So basically implementing the equations of motion (the Kuramoto dynamics equation) for both oscillators.

In general, the equations of motion are:

$$\frac{\mathrm{d}\theta_i}{\mathrm{d}t} = \omega_i + K \sum_{j}^{N}A_{ij}\sin(\theta_j-\theta_i)$$

then for two oscillators $\theta_1$ and $\theta_2$ we have:

$$\frac{\mathrm{d}\theta_1}{\mathrm{d}t} = \omega_1 + K\sin(\theta_2-\theta_1)$$
$$\frac{\mathrm{d}\theta_2}{\mathrm{d}t} = \omega_2 + K\sin(\theta_1-\theta_2)$$

You can start by setting all the parameters to be simple numbers, for instance two recognizable natural frequencies, let's say $\omega_1 = 1$, $\omega_2 = 4$, such that one is faster than the other. And no interaction, $i.e.,$ coupling $K = 0$.
 
Before doing any analysis let's see how the phases behave without interactions:
"""

# ╔═╡ 04c36bfd-5aee-42e2-b110-9e376881a91c
md"""
Select $\omega_{1}$: $(@bind ω_1_tutorial_1 confirm(NumberField(0.00:0.1:10, default=1.0)))

Select $\omega_{2}$: $(@bind ω_2_tutorial_1 confirm(NumberField(0.00:0.1:10, default=4.0)))
"""

# ╔═╡ 643768b0-234d-11ef-1090-bd0594a6fa7e
let
	A = [0; 1;; 1; 0]
	N = size(A)[1]
	Δt = 1e-3
	sim_time = 20
	seedval = 123
	ω = [ω_1_tutorial_1,ω_2_tutorial_1]
	K = 0
	A = A .* K

	yticks = ([0, π, 2*π], [L"0", L"\pi", L"2\pi"])

	xticks = ([0, π/2, π, 3*π/2, 2*π, 3*π, 4*π, 5*π, 6*π], [L"0", L"\frac{\pi}{2}", L"\pi", L"3\frac{\pi}{2}", L"2\pi", L"3\pi", L"4\pi", L"5\pi", L"6\pi"])
	
	# simulations settings
	Random.seed!(seedval)
	
	# simulation parameters
	steps = (0.0+Δt):Δt:sim_time
	no_steps = length(steps)
	
	# storing parameters
	save_ratio = 1
	no_saves = round(Integer, no_steps / save_ratio)
	θ_now = [0,0]  # random init conditions
	θs = zeros(no_saves, N)
	θs[1, :] = θ_now
	
	save_counter = 1
	
	for t in 2:no_steps
	    
	    # update phases
	    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)')
	    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 
	
	    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
	    θ_now += Δt .* (ω + k1)
	    save_counter += 1
	
	    # save θ
	    if save_counter % save_ratio == 0
	        θs[round(Integer, save_counter / save_ratio), :] = θ_now
	    end
	    
	end

	plot(macro_op(θs))

	plt1 = plot(steps, θs[:, 1] .% (2π), label = L"i = 1",
		ylabel = L"\theta_{i}",
		xticks = xticks,
		yticks = yticks,
	)
	plot!(steps, θs[:, 2] .% (2π), label = L"i = 2")

	plt2 = plot(steps, sin.(θs[:, 1]),
		label = L"i = 1",
		ylabel = L"\theta_{i}", xticks = xticks
	)
	plot!(steps, sin.(θs[:, 2]), label = L"i = 2", yticks = ([-1,0,1], [L"-1", L"0", L"1"]))

	plot(
		[plt1, plt2]...,
		layout = (2,1),
		size = (700, 350),
		tickfontsize = 12,
		legendfontsize = 8,
		legend = :topright,
		xlabelfontsize = 12,
		ylabelfontsize = 12,
		left_margin = 8Plots.mm,
		right_margin = 8Plots.mm,
		bottom_margin = 8Plots.mm,
		top_margin = 8Plots.mm
	)
end

# ╔═╡ 47159082-f661-4351-bcf9-fc64fa16a755
md"""
here I plotted the raw phases data above (with modulus $2\pi$), and below I plotted the sine of of the phases. So if the phases are $\theta_1(t)$ and $\theta_2(t)$, I plotted the $\sin(\theta_1(t))$.
 
Then try setting $K>0$ and compare the phases plots:

Here I used $K = 1$ (you can change this using the button below) and you can see that the fast oscillator's dynamic doesn't change much, while the slower oscillator's does.
"""

# ╔═╡ 17d1a4a6-7545-4a6f-8eea-b943322920ce
md"""
Select $K$: $(@bind K_tutorial_1 confirm(NumberField(0.00:0.1:4, default=1.0)))
"""

# ╔═╡ eff4e3e4-f6a8-46a3-8c19-c811d97301b4
let
	A = [0; 1;; 1; 0]
	N = size(A)[1]
	Δt = 1e-3
	sim_time = 10
	seedval = 123
	ω = [1,4]
	K = K_tutorial_1
	A = A .* K
	
	# simulations settings
	Random.seed!(seedval)
	
	# simulation parameters
	steps = (0.0+Δt):Δt:sim_time
	no_steps = length(steps)
	
	# storing parameters
	save_ratio = 1
	no_saves = round(Integer, no_steps / save_ratio)
	θ_now = [0,0]  # random init conditions
	θs = zeros(no_saves, N)
	θs[1, :] = θ_now
	
	save_counter = 1
	
	for t in 2:no_steps
	    
	    # update phases
	    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)')
	    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 
	
	    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
	    θ_now += Δt .* (ω + k1)
	    save_counter += 1
	
	    # save θ
	    if save_counter % save_ratio == 0
	        θs[round(Integer, save_counter / save_ratio), :] = θ_now
	    end
	    
	end

	plot(macro_op(θs))

	plt1 = plot(steps, θs[:, 1] .% (2π), label = "oscillator 1")
	plot!(steps, θs[:, 2] .% (2π), label = "oscillator 2")

	plt2 = plot(steps, sin.(θs[:, 1]), label = "oscillator 1")
	plot!(steps, sin.(θs[:, 2]), label = "oscillator 2")

	plot(
		[plt1, plt2]...,
		layout = (2,1),
		size = (700, 350)
	)
end

# ╔═╡ f4af7d76-d81e-4aaa-808f-0a1c0089b238
md"""
Very cool everything seems to work nicely, now let's try with a few more oscillators
"""

# ╔═╡ 8f38dded-217b-4d5d-adb4-819d9d7eaee9
md"""
### Experiment with large network
"""

# ╔═╡ 250df70c-7156-44b9-92a9-df1f976acfbf
md"""
# Functions
"""

# ╔═╡ eed3e5e6-2b1b-465e-8995-82733100483a
"""
	function _partition_mat(n, B)

	Inputs:
	
		n::AbstractArray 		vector of length 3, number of l-1 blocks in each l block
		B::AbstractArray 		number of blocks in layer l

	Output:

		P::AbstractMatrix 		partition matrix

"""
function _partition_mat(n, B)

	N = prod(n)
	L = length(n)
	P = zeros(L-1, N)
	
	for i in 1:(L-1)
		P[i, :] = vcat(
			[repeat([j], prod(n[1:i])) for j in 1:B[i]]...
		)
	end

	return P
	
end

# ╔═╡ f18439fb-08b0-49ef-9676-e5a436611f02
let
	# simulations settings
	seedval = 123
	Random.seed!(seedval)
	
	# number of blocks in each layer
	B = [32, 8, 1]
	# number of oscillators per module
	n = [16, 4, 8]
	
	# construct 3-layer network
	N = prod(n)
	A = zeros(N, N)
	P = _partition_mat(n, B)
	p = [0.95, 0.4, 0.15] # just a test
	for i in 1:N
		for j in i+1:N
			if P[1, i] == P[1, j] && rand() < p[1]
				A[i,j] = 1
				A[j,i] = 1
			elseif P[2, i] == P[2, j] && rand() < p[2]
				A[i,j] = 1
				A[j,i] = 1
			elseif rand() < p[3]
				A[i,j] = 1
				A[j,i] = 1
			end
		end
	end

	println(mean([count(x -> x==1, A[i, :]) for i in 1:size(A)[1]]))
	
	# matrix plot
	plt1 = heatmap(A, yflip = true)

end

# ╔═╡ c04b9210-4fb1-4dde-85c6-8bacdd21967c
"""
	function Laplacian(A::AbstractMatrix; method = :graph)

Function used to compute the Laplacian of the graph defined by the
adjacency matrix `A`.

If `method == :graph` --> compute the graph Laplacian;

If `method == :normalized` --> compute the normalized graph Laplacian.

"""
function Laplacian(A::AbstractMatrix; method = :graph)

	N = length(A[1, :])
	
	if method == :graph
		
		D = diagm([sum(A[i, :]) for i in 1:N])
		
		return D - A
		
	elseif method == :normalized

		𝐼 = diagm(ones(N))
		k = [sum(A[i, :]) for i in 1:N]
		A′ = zeros(N, N)
		
		for i in 1:N
			for j in 1:N
				if A[i,j] != 0
					A′[i,j] = A[i,j] / √(k[i]*k[j])
				end
			end
		end

		return 𝐼 - A′
	end
	
end

# ╔═╡ 47f97a83-8490-4e49-a3f7-3e239862ce9f
let
	# simulations settings
	seedval = 123
	Random.seed!(seedval)
	
	# number of blocks in each layer
	B = [32, 8, 1]
	# number of oscillators per module
	n = [16, 4, 8]
	
	# construct 3-layer network
	N = prod(n)
	A = zeros(N, N)
	P = _partition_mat(n, B)
	p = [0.95, 0.4, 0.15] # just a test
	for i in 1:N
		for j in i+1:N
			if P[1, i] == P[1, j] && rand() < p[1]
				A[i,j] = 1
				A[j,i] = 1
			elseif P[2, i] == P[2, j] && rand() < p[2]
				A[i,j] = 1
				A[j,i] = 1
			elseif rand() < p[3]
				A[i,j] = 1
				A[j,i] = 1
			end
		end
	end

	# matrix plot
	plt1 = heatmap(A, yflip = true)

	# get Laplacian
	L = Laplacian(A)
	λs, v = eigen(L)

	# plot Laplacian eigenvalues
	plt2 = plot(
		1 ./ λs[2:end],
		collect(2:lastindex(λs)),
		xaxis = :log10,yaxis = :log10,
		marker = true,
		markersize = 5
	)

	# simulation parameters
	Δt = 1e-3
	sim_time = 10
	ω = repeat([1], N)
	K = 1
	A = A .* K
	
	# simulation parameters
	steps = (0.0+Δt):Δt:sim_time
	no_steps = length(steps)
	
	# storing parameters
	save_ratio = 1
	no_saves = round(Integer, no_steps / save_ratio)
	θ_now = rand(Uniform(-π, π), N)  # random init conditions
	θs = zeros(no_saves, N)
	θs[1, :] = θ_now
	α = π/2 - 0.1
	
	save_counter = 1
	
	for t in 2:no_steps
	    
	    # update phases
	    θj_θi_mat = (repeat(θ_now',N) - repeat(θ_now',N)') .- α
	    setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero 
	
	    k1 = map(sum, eachrow(A .* sin.(θj_θi_mat)))
	    θ_now += Δt .* (ω + k1)
	    save_counter += 1
	
	    # save θ
	    if save_counter % save_ratio == 0
	        θs[round(Integer, save_counter / save_ratio), :] = θ_now
	    end
	    
	end

	plt3 = plot(macro_op(θs))

	global θs_test = θs
end

# ╔═╡ 3bf4e87c-5702-4948-b03b-c6be7340c55a
let
	# simulations settings
	seedval = 123
	Random.seed!(seedval)
	
	# number of blocks in each layer
	B = [24, 8, 1]
	# number of oscillators per module
	n = [16, 3, 8]
	
	# construct 3-layer network
	N = prod(n)
	A = zeros(N, N)
	P = _partition_mat(n, B)
	p = [0.9, 0.4, 0.1] # just a test
	
	θs = θs_test

	# simulation parameters
	Δt = 1e-3
	sim_time = 10
	ω = repeat([1], N)
	K = 1
	A = A .* K
	
	# simulation parameters
	steps = (0.0+Δt):Δt:sim_time
	no_steps = length(steps)

	pops_KOP = zeros(length(steps), B[2])
	for i in 1:B[2]
		pops_KOP[:, i] = macro_op(θs[:, n[1]*(i-1)+1:n[1]*i])
	end

	plot(pops_KOP[end-3000:end, 8], lw = 5, ylims = (0,1))
	# plot(pops_KOP[end-2000:end-1000, :], lw = 5)

	# γ = 0.8
	# pops_coals = zeros(size(pops_KOP))
	# for i in 1:length(steps)
	# 	for j in 1:B[2]
	# 		if pops_KOP[i,j] > γ
	# 			pops_coals[i,j] = 1
	# 		end
	# 	end
	# end

	# heatmap(pops_coals)
	# pops_coals
end

# ╔═╡ Cell order:
# ╟─e5277ce7-db25-4a0e-b2fe-f150a0a5ca3b
# ╟─15f1cffb-e9b8-4abb-892e-48b4e7bea33c
# ╟─017e5156-b4f2-4c90-9bf4-352c7d877ed1
# ╟─0f884c5e-a00e-4dda-95ce-89d3d85cf94c
# ╟─04c36bfd-5aee-42e2-b110-9e376881a91c
# ╟─643768b0-234d-11ef-1090-bd0594a6fa7e
# ╟─47159082-f661-4351-bcf9-fc64fa16a755
# ╟─17d1a4a6-7545-4a6f-8eea-b943322920ce
# ╟─eff4e3e4-f6a8-46a3-8c19-c811d97301b4
# ╟─f4af7d76-d81e-4aaa-808f-0a1c0089b238
# ╟─8f38dded-217b-4d5d-adb4-819d9d7eaee9
# ╟─47f97a83-8490-4e49-a3f7-3e239862ce9f
# ╟─f18439fb-08b0-49ef-9676-e5a436611f02
# ╟─3bf4e87c-5702-4948-b03b-c6be7340c55a
# ╟─250df70c-7156-44b9-92a9-df1f976acfbf
# ╟─eed3e5e6-2b1b-465e-8995-82733100483a
# ╟─c04b9210-4fb1-4dde-85c6-8bacdd21967c

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

# ╔═╡ 0f884c5e-a00e-4dda-95ce-89d3d85cf94c
md"""
# Tutorial 1

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


# ╔═╡ Cell order:
# ╟─15f1cffb-e9b8-4abb-892e-48b4e7bea33c
# ╟─017e5156-b4f2-4c90-9bf4-352c7d877ed1
# ╟─0f884c5e-a00e-4dda-95ce-89d3d85cf94c
# ╟─04c36bfd-5aee-42e2-b110-9e376881a91c
# ╟─643768b0-234d-11ef-1090-bd0594a6fa7e
# ╟─47159082-f661-4351-bcf9-fc64fa16a755
# ╟─17d1a4a6-7545-4a6f-8eea-b943322920ce
# ╟─eff4e3e4-f6a8-46a3-8c19-c811d97301b4
# ╟─f4af7d76-d81e-4aaa-808f-0a1c0089b238
# ╠═8f38dded-217b-4d5d-adb4-819d9d7eaee9

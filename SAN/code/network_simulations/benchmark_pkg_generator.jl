import Pkg
Pkg.generate("joint_timing")
Pkg.activate("joint_timing")
Pkg.add("Test")
Pkg.add("BenchmarkTools")
Pkg.add("Profile")
Pkg.add(Pkg.PackageSpec(;name="Distributions", version="0.24.15"))
Pkg.add("Random")
Pkg.add("SpecialFunctions")
Pkg.add("LinearAlgebra")
Pkg.add("ForwardDiff")
Pkg.add("FiniteDiff")
Pkg.add("Zygote")
Pkg.add("DistributionsAD")
Pkg.add("JuMP")
Pkg.add("Ipopt")
Pkg.add("NLsolve")
Pkg.add("DataFrames")
Pkg.add("XLSX")
Pkg.add("Missings")
Pkg.add("JLD2")
Pkg.add("CSV")
Pkg.add("Quadrature")
Pkg.add("Cuba")
Pkg.add("Cubature")
Pkg.add("HCubature")

using Test, BenchmarkTools, Profile
using Random, SpecialFunctions, LinearAlgebra, ForwardDiff, FiniteDiff, Zygote, Distributions, DistributionsAD
using JuMP, Ipopt, NLsolve
using DataFrames, XLSX, Missings, JLD2, CSV
using Quadrature, Cuba, Cubature, HCubature
using CSV, DataFrames
using CairoMakie, LaTeXStrings
using Statistics
using LsqFit

colors = [:red, :blue, :green, :orange, :purple, :brown, :pink, :gray, :black]
# colors = ["#4477AA", "#E36D44", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
markersize = 12
linestyle = [:solid, :dash, :dot, :dashdot, :dashdotdot]
linewidth = 2
markerstyle = [:rect, :circle, :diamond, :utriangle, :dtriangle, :rtriangle, :ltriangle]
strokecolor = :black
strokewidth = 1.0

function geometric_mean(x)
    return prod(Float64.(x))^(1/length(x))
end

function max_n(x, n)
    # get the n largest elements in x
    xs = sort(x)
    return xs[end-n+1:end]
end

hstyle = :dot
hwidth = 1.5

function geometric_mean(x)
    return prod(Float64.(x))^(1/length(x))
end
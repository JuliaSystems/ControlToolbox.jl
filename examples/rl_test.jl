using ControlCore
using ControlToolbox

using PyPlot

"""
    `plot_rlocus(sys, [K])`
    Plots the root locus of `sys`.

    If `K` is omitted, the pole trajectories are plotted as evenly spaced
    samples.
"""
function plot_rlocus(sys::ControlCore.LtiSystem)
    plist = rlocus(G)[1]

    fig = figure() #"rlocus_plot", figsize=(10,10))
    ax = axes()

    scatter(real(plist),imag(plist),marker="o")
    title("Root Locus")
    xlabel("Re")
    ylabel("Im")
    grid("on")
end


G = tf([2.0, 5, 1],[1.0, 2, 3])

rlocus(G)
rlocus(G, 3)[1]
rlocus(G, [1, 2, 5, 100])

plot_rlocus(G)

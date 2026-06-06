#=

This file makes plots, tidying up after arxiv v1
April-May 2026

=#

using ToSteerOrNot, Plots, OneTwoFive, Statistics, Random, DelimitedFiles, Pkg
SAVE = false


####################
##### FIGURE 1 #####


matt_g_i_v = DelimitedFiles.readdlm(joinpath(Pkg.pkgdir(ToSteerOrNot), "docs/Analytic_Frontier_Points.csv"), ' ')

function myscatter!(xs::Vector; m=:circle, c=colour2(xs[1]), lab="", kw...)
    scatter!(getfield.(xs, :info), getfield.(xs, :speed); m, c, lab, kw...)
end
function myscatter!(p::Plots.Plot, xs::Vector; m=:circle, c=colour2(xs[1]), lab="", kw...)
    scatter!(p, getfield.(xs, :info), getfield.(xs, :speed); m, c, lab, kw...)
end
function myplot!(xs::Vector; c=colour2(xs[1]), lab="", kw...)
    plot!(getfield.(xs, :info), getfield.(xs, :speed); c, lab, kw...)
end


tumble40 = [max2Wcostumble(γ) for γ in logrange(0.001, 2, 40)]

tumble20strong = [max2Wstrong(@show(γ), 180) for γ in logrange(0.0001, 0.005, 30)];


reverse30 = [max2Wcosreverse(γ) for γ in logrange(0.001, 1.5, 30)]
reverse20strong = [max2Wstrongreverse(γ, 180) for γ in logrange(0.0001, 0.01, 20)]

unsigned40 = [max2Wcosabs(γ) for γ in logrange(0.0005, 2, 40)]
unsigned50 = [max2Wcosabs(γ) for γ in logrange(0.0002, 2, 50)]


#===== Talk =====#

fig1_one = let
    plot(xguide="information rate i/Dᵣ", yguide="climbing speed v/v₀")

    # tumble
    c=COLOURS2D[:lambda]
    plot!(x -> sqrt(x/4), 1e-4, 3; lab="", s=:dash, c)
    annotate!(0.1, 0.2, Plots.text("sqrt(i/4Dᵣ)", 10, c, rotation=45))

    myscatter!(tumble40, lab="tumble λ(θ)")
    # myscatter!(tumble20strong, m=:square, alpha=0.5, ms=3)

    # axes
    plot!(xaxis=:log10, yaxis=:log10, size=(400, 280))
    plot!(legend=false)
    # plot!(legendfontsize=10, legend=:bottomright, foreground_color_legend=nothing)
    plot!(xlims=[0.03, 50], xitcks=[0.1, 1, 10, 100], ylims=[0.1, 1.05], yticks=[0.1, 1])
    # plot!(xitcks=[0.1, 1, 10, 100])  # hmm
    plot!(dpi=300)

end

SAVE && savefig(fig1_one, "april01-fig1one.png")


fig1_two = let
    plot(xguide="information rate i/Dᵣ", yguide="climbing speed v/v₀")

    # steering
    plot!(matt_g_i_v[:,3], matt_g_i_v[:,2], c=:red, lab="", l=2, fill=1, fillalpha=0.1, fillcolor=:grey)
    plot!([NaN, NaN], c=:red, l=2, lab="steer μ(θ)")

    plot!(x -> sqrt(x/2), 1e-4, 1; lab="", s=:dash, c=:red)
    annotate!(0.07, 0.21, Plots.text("sqrt(i/2Dᵣ)", 10, :red, rotation=45))

    annotate!(0.2, 0.6, Plots.text("impossible", 10, :grey, rotation=45)) # "impossible:\nv > vₛₜₑₑᵣ(i)"

    # tumble
    c=COLOURS2D[:lambda]
    plot!(x -> sqrt(x/4), 1e-4, 1; lab="", s=:dash, c)
    annotate!(0.1, 0.18, Plots.text("sqrt(i/4Dᵣ)", 10, c, rotation=45))

    myscatter!(tumble40, lab="tumble λ(θ)")
    # myscatter!(tumble20strong, m=:square, alpha=0.5, ms=3)

    # axes
    plot!(xaxis=:log10, yaxis=:log10, size=(400, 280))
    plot!(legendfontsize=10, legend=:bottomright, foreground_color_legend=nothing)
    plot!(xlims=[0.03, 50], xitcks=[0.1, 1, 10, 100], ylims=[0.1, 1.05], yticks=[0.1, 1])
    plot!(xitcks=[0.1, 1, 10, 100])  # hmm
    plot!(dpi=300)

end

SAVE && savefig(fig1_two, "april01-fig1two.png")


fig1_three = let p = deepcopy(fig1_two)
    # unsigned
    c = COLOURS2D[:abs]
    plot!(p, x -> sqrt(x/8), 1e-4, 1; lab="", s=:dash, c)
    annotate!(p, 0.15, 0.155, Plots.text("sqrt(i/8Dᵣ)", 9, c, rotation=45))

    # myscatter!(unsigned40, m=:downtri, lab="steer μₙₒ(|θ|)")
    myscatter!(unsigned40, m=:downtri, lab="steer μ(θ)=μ(-θ)")
end

SAVE && savefig(fig1_three, "april01-fig1three.png")



fig1_rev = let p = deepcopy(fig1_two)
    # reverse
    c = COLOURS2D[:zeta]
    plot!(p, x -> 2/pi, 3, 1000; lab="", s=:dash, c)
    annotate!(p, 30, 0.59, Plots.text("2/π", 9, c))

    myscatter!(p, reverse20strong, m=:square, alpha=0.5, ms=3)
    myscatter!(p, reverse30, m=:uptri, lab="reverse λᵣₑᵥ(θ)", ms=5)

end

SAVE && savefig(fig1_rev, "april01-fig1rev.png")




#===== Paper =====#

# These are for Fig3, moved up because I added them to Fig1 to see...
pair1 = max2Wpair(0.5, 36; log=true)
pair2 = max2Wpair(0.05, 72; log=true)
pair3 = max2Wpair(0.005, 72; log=true, iterations=3000, g_abstol=1e-10)

# These are
_jump125log = @time [try
        max2Wjump(@show(γ); mirror=false, log=true)
    catch e
        @show nothing
    end for γ in onetwofive(0.001, 10.0)];
jump125log = filter(!isnothing, _jump125log)
whole2heatmap(jump125log, title="any jump, working in log")



fig1 = let
    plot(xguide="information rate i/Dᵣ", yguide="climbing speed v/v₀")

    # annotate!(0.2, 0.7, Plots.text("impossible", 9, :grey, rotation=45)) # "impossible:\nv > vₛₜₑₑᵣ(i)"
    annotate!(0.2, 0.7, Plots.text("forbidden", 9, :grey, rotation=45))

    # steering
    plot!(matt_g_i_v[:,3], matt_g_i_v[:,2], c=:red, lab="", l=2, fill=1, fillalpha=0.1, fillcolor=:grey)
    plot!([NaN, NaN], c=:red, l=2, lab="steer μ(θ), odd, exact soln.")

    # plot!(x -> 1-1/(2x), 0.5, 1000; lab="", s=:dot, c=:red)
    # annotate!(0.5, 0.15, Plots.text("1-Dᵣ/2i", 9, :red, rotation=82))

    # reverse
    c = COLOURS2D[:zeta]
    plot!(x -> sqrt(x/2), 1e-4, 1; lab="", s=:dash, c)
    plot!(x -> 2/pi, 3, 1000; lab="", s=:dash, c)
    annotate!(0.07, 0.21, Plots.text("sqrt(i/2Dᵣ)", 9, c, rotation=45))
    # annotate!(0.07, 2/pi, Plots.text("2/π", 9, c))
    annotate!(30, 0.59, Plots.text("2/π", 9, c))

    myscatter!(reverse30, m=:uptri, lab="reverse λᵣₑᵥ(θ), even")
    myscatter!(reverse20strong, m=:square, alpha=0.5, ms=3)

    # tumble
    c=COLOURS2D[:lambda]
    plot!(x -> sqrt(x/4), 1e-4, 1; lab="", s=:dash, c)
    annotate!(0.1, 0.18, Plots.text("sqrt(i/4Dᵣ)", 9, c, rotation=45))

    myscatter!(tumble40, lab="tumble λ(θ), even") # , delay τ=0")
    myscatter!(tumble20strong, m=:square, alpha=0.5, ms=3)

    # unsigned
    c = COLOURS2D[:abs]
    plot!(x -> sqrt(x/8), 1e-4, 1; lab="", s=:dash, c)
    annotate!(0.15, 0.155, Plots.text("sqrt(i/8Dᵣ)", 9, c, rotation=45))

    myscatter!(unsigned40, m=:downtri, lab="steer μ(|θ|) constrained even")

    # jump
    myscatter!(jump125log, lab="any λ(Δθ,θ), ≈ steering", m=:diamond, ms=3, c=:red)

    myscatter!([pair1, pair2, pair3], c=:black, m=:star, lab="constrained λ(|Δθ|,θ), even", ms=5)

    # axes
    plot!(xaxis=:log10, yaxis=:log10)
    plot!(xlims=[0.03, 50], ylims=[0.1, 1.05])
    # plot!(yticks=([0.1, 0.3, 1], ["0.1", "0.3", "1"]))
    # plot!(yticks=([0.1, 0.2, 0.5, 1], ["0.1", "0.2", "0.5", "1"]))
    plot!(xitcks=([0.1, 1, 10], ["0.1", "1", "10"]))  # ignored?
    plot!(legendfontsize=9)
    plot!(minorgrid=true)
    plot!(size=(550, 380))
end

SAVE && savefig(fig1, "april01-fig1.pdf")



####################
##### FIGURE 2 #####

# Same solutions as in nov-01

tum001 = max2Wcostumble(2.47)
rev001 = max2Wcosreverse(3.5)

ste001 = max2Wsinsteer(3.45)
abs001 = max2Wcosabs(1.75)


tum01 = max2Wcostumble(0.76)
rev01 = max2Wcosreverse(1.03, 36)

ste01 = max2Wsinsteer(1.03, 36)
# abs01 = max2Wcosabs(0.53)
abs01 = max2Wcosabs(0.55)

tum_strong = max2Wstrong(0.001, 144)
rev_strong = max2Wstrongreverse(0.00058, 150)

ste_fast = max2Wsinsteer(0.15, 36)
steertojump(ste_fast, 2) |> whole2heatmap
abs_fast = max2Wcosabs(0.05)


#===== Talk =====#

fig2A_heat = let
    whole2heatmap(tum01; color=cgrad([:lightgrey, :lightgrey, :blue, :darkblue]), sigdigits=2, both=false)
    plot!(xguide="heading θ", yguide="turn angle Δθ", colorbar_title="\nrate λ(Δθ,θ)")
    plot!(clims=(0, 0.028))
    plot!(size=(300,250), dpi=300)
    plot!(margin=0*Plots.mm, right_margin=4*Plots.mm)
end

SAVE && savefig(fig2A_heat, "april01-fig2Aheat.png")

fig2B_heat = let
    whole2heatmap(rev01; color=cgrad([:lightgrey, :green, :darkgreen]), sigdigits=2, both=false)
    plot!(xguide="heading θ", yguide="turn angle Δθ", colorbar_title="\nrate λ(Δθ,θ)")
    plot!(size=(300,250), dpi=300)
    plot!(margin=0*Plots.mm, right_margin=4*Plots.mm)
end

SAVE && savefig(fig2B_heat, "april01-fig2Bheat.png")

fig2AB_heat = let
    p1 = deepcopy(fig2A_heat)
    plot!(p1, colorbar_title="")
    p2 = deepcopy(fig2B_heat)
    plot!(p2, xguide="θ", yguide="", colorbar_title="rate λ(Δθ,θ)")
    plot(p1, p2, size=(500,220))
    plot!(margin=0*Plots.mm, left_margin=2*Plots.mm, right_margin=0*Plots.mm, bottom_margin=3*Plots.mm)
end

SAVE && savefig(fig2AB_heat, "april01-fig2ABheat.png")


# let
#     whole2heatmap(ste01; color=cgrad([:lightgrey, :green, :darkgreen]), sigdigits=2, both=false)
#     plot!(xguide="heading θ", yguide="turn angle Δθ", colorbar_title="rate λ(Δθ,θ)")
# end



fig2Avert = let
    p1 = whole2plot(tum001; sigdigits=2, yguide="", yguide2="tumble rate λ(θ)")
    plot!(p1, xguide="")
    # plot!(p1, left_margin=0*Plots.mm)

    p2 = whole2plot(tum_strong; sigdigits=2, yguide="steady-state p(θ)", xguide="heading θ", yticks=[0, 1000])
    plot!(p2, yguidefontcolor=:darkgrey)
    # plot!(p2, right_margin=0*Plots.mm, left_margin=1*Plots.mm)


    plot(p1, p2,
        layout=(2,1),
        size=(250,500),
        # left_margin=Plots.mm*4, right_margin=Plots.mm*5,  # unforunately apply to all!
        # top_margin=10*Plots.mm,
        plot_title=" ",
        margin=Plots.mm*0, # bottom_margin=Plots.mm*0,
        dpi=300,
    )
end

SAVE && savefig(fig2Avert, "april01-fig2Avert.png")

fig2A = let
    p1 = whole2plot(tum001; sigdigits=2, xguide="heading θ", yguide="steady-state p(θ)")
    plot!(p1, yguidefontcolor=:darkgrey)
    plot!(p1, left_margin=4*Plots.mm)

    p2 = whole2plot(tum_strong; sigdigits=2, yguide2="tumble rate λ(θ)", xguide="θ", yticks=[0, 1000])
    # plot!(p2, right_margin=0*Plots.mm, left_margin=1*Plots.mm)

    plot(p1, p2,
        layout=(1,2),
        size=(500,250),
        dpi=300,
    )
end

SAVE && savefig(fig2A, "april01-fig2A.png")

fig2B = let
    p4 = whole2plot(rev001, sigdigits=2, xguide="heading θ", yguide="steady-state p(θ)", yticks=[0, 0.25, 0.5])
    plot!(p4, yguidefontcolor=:darkgrey)
    plot!(p4, left_margin=5*Plots.mm)

    p5 = whole2plot(rev_strong, sigdigits=2, xguide="θ", yguide2="reverse rate λᵣₑᵥ(θ)", yticks=[0, 200, 400])
    plot!(p5, right_margin=2*Plots.mm, left_margin=0*Plots.mm)

    plot(p4, p5,
        layout=(1,2),
        size=(500,250),
        dpi=300,
    )
end

SAVE && savefig(fig2B, "april01-fig2B.png")

p7 = whole2plot(ste001; sigdigits=2, yguide2="steering rate μ(θ)", xguide="heading θ")
plot!(p7, dpi=300, size=(250, 250), right_margin=0*Plots.mm, bottom_margin=0*Plots.mm)

SAVE && savefig(p7, "april01-fig2Cone.png")

fig2Ctwin = let sol = ste_fast
    p7 = whole2plot(sol; sigdigits=2, xguide="heading θ", yguide="steady-state p(θ)")
    plot!(p7, yguidefontcolor=:darkgrey, left_margin=4*Plots.mm, right_margin=-2*Plots.mm, title="")

    p8 = whole2heatmap(steertojump(sol, 2); sigdigits=2, both=false, color=color=cgrad([:lightgrey, :red, :darkred]))
    # _cmax = 100 # maximum(tmp.lambda)
    plot!(p8, xguide="θ", yguide="turn angle Δθ", #colorbar_title="steering rate μ(θ)", colorbar_titlefontcolor=:red,
        colorbar_ticks=[0,1], xlims=[0,2pi],)
    # plot!(twinx(p8), [NaN], [NaN]; lab="", xlims=[0,2pi], yticks=[], yguide="steering rate μ(θ)/Dᵣ\n(using sign)", yguidefontcolor=:red)
    plot!(p8; title="")

    v = round(sol.speed; sigdigits=2)
    i = round(sol.info; sigdigits=2)
    plot_title = "i/Dᵣ=$i, v/v₀=$v, and α=20°"

    plot(p7, p8;
        layout=(1,2), plot_title,
        size=(500,250),
        dpi=300,
    )
end

SAVE && savefig(fig2Ctwin, "april01-fig2Ctwin.png")

fig2D = let
    p9 = whole2plot(abs001, sigdigits=2, yguide="steady-state p(θ)", xguide="heading θ")
    plot!(p9, yguidefontcolor=:darkgrey)
    plot!(p9, left_margin=5*Plots.mm, right_margin=0*Plots.mm)

    p10 = whole2plot(abs_fast, sigdigits=2, xguide="θ", yguide2="steering rate μ(θ)=μ(-θ)")

    # p4 = whole2plot(rev001, sigdigits=2, xguide="heading θ", yguide="p(θ)", yticks=[0, 0.25, 0.5])
    # plot!(p4, yguidefontcolor=:darkgrey)
    # plot!(p4, left_margin=5*Plots.mm)

    # p5 = whole2plot(rev_strong, sigdigits=2, xguide="θ", yguide2="reverse rate λ(θ)", yticks=[0, 200, 400])
    # plot!(p5, right_margin=2*Plots.mm, left_margin=0*Plots.mm)

    plot(p9, p10,
        layout=(1,2),
        size=(500,250),
        dpi=300,
    )
end

SAVE && savefig(fig2D, "april01-fig2D.png")


#===== Paper =====#

fig2AB = let
    # First row
    p1 = whole2plot(tum001; sigdigits=2, yguide="steady-state p(θ)")
    plot!(p1, xguide="", yguidefontcolor=:darkgrey)
    plot!(p1, left_margin=1*Plots.mm)

    p2 = whole2plot(tum_strong; sigdigits=2, xguide="", yticks=[0, 500, 1000], yguide2="tumble rate λ(θ)/Dᵣ")
    plot!(p2, right_margin=0*Plots.mm, left_margin=1*Plots.mm)

    # Second row
    p4 = whole2plot(rev001, sigdigits=2, xguide="heading θ", yguide="p(θ)", yticks=[0, 0.25, 0.5])
    plot!(p4, yguidefontcolor=:darkgrey)
    plot!(p4, left_margin=1*Plots.mm)

    p5 = whole2plot(rev_strong, sigdigits=2, xguide="θ", yticks=[0, 200, 400], yguide2="reverse rate λᵣₑᵥ(θ)/Dᵣ")
    plot!(p5, right_margin=0*Plots.mm, left_margin=1*Plots.mm)

    plot(p1, p2, p4, p5;
        layout=(2,2),
        size=(500,500),
        bottom_margin=Plots.mm*0,
    )
end


fig2CD = let
    # First row
    p7 = whole2plot(ste001; sigdigits=2, yguide=" ", yguide2="steering rate μ(θ)/Dᵣ")
    plot!(p7, xguide="", right_margin=-5*Plots.mm)
    plot!(p7, left_margin=1*Plots.mm)

    p8 = plot(axis=false, ticks=[])

    # Second row
    p9 = whole2plot(abs001, sigdigits=2, xguide="θ", yguide=" ")
    plot!(p9, yguidefontcolor=:darkgrey, right_margin=-5*Plots.mm)
    plot!(p9, left_margin=1*Plots.mm)

    p10 = whole2plot(abs_fast, sigdigits=2, xguide="θ", yguide2="steering rate μ(|θ|)/Dᵣ")
    plot!(p10, right_margin=0*Plots.mm, left_margin=1*Plots.mm)

    plot(p7, p8,
        p9, p10,
        # colorbar=false,
        layout=(2,2),
        size=(480, 500),
        bottom_margin=Plots.mm*0,
    )
end


fig2Dheat = let

    p8 = whole2heatmap(steertojump(ste_fast, 2); sigdigits=2, both=false, color=color=cgrad([:lightgrey, :red, :darkred]))
    # _cmax = 100 # maximum(tmp.lambda)
    plot!(p8, xguide="", yguide="turn angle Δθ", colorbar_title="\nλ(Δθ,θ)/Dᵣ", colorbar_titlefontcolor=:black, colorbar_ticks=[0,1], xlims=[0,2pi],)
    # plot!(twinx(p8), [NaN], [NaN]; lab="", xlims=[0,2pi], yticks=[], yguide="steering rate μ(θ)/Dᵣ\n(using sign)", yguidefontcolor=:red)
    plot!(p8, right_margin=5*Plots.mm, xguide=" ", left_margin=-2*Plots.mm)

    plot!(p8, size=(250, 250), bottom_margin=Plots.mm*0)
end

SAVE && savefig(fig2AB, "april01-fig2AB.svg")
SAVE && savefig(fig2CD, "april01-fig2CD.svg")
SAVE && savefig(fig2Dheat, "april01-fig2Dheat.svg")


###################################
##### FIGURE 2 - TRAJECTORIES #####

#===== Talk =====#

let sol = tum01
    Random.seed!(1); trajplot(sol, 0:0.01:50, 1; limit=5, lab="", ms=5, msw=2, l=2)
    plot!(xguide="", yguide="", margin=0*Plots.mm, aspect_ratio=1)
    plot!(xticks=[], yticks=[])
    annotate!(0, 0, Plots.text("  start", 11, :black, :left))
    v = round(sol.speed; sigdigits=2)
    i = round(sol.info; sigdigits=2)
    plot!(title = "i/Dᵣ=$i, v/v₀=$v", dpi=300)
end

SAVE && savefig("april01-trajA.png")

let sol = rev01
    Random.seed!(2); trajplot(sol, 0:0.01:50, 1; limit=5, lab="", ms=5, msw=2, l=2)
    plot!(xguide="", yguide="", margin=0*Plots.mm, aspect_ratio=1)
    plot!(xticks=[], yticks=[])
    annotate!(0, 0, Plots.text("  start", 11, :black, :left))
    v = round(sol.speed; sigdigits=2)
    i = round(sol.info; sigdigits=2)
    plot!(title = "i/Dᵣ=$i, v/v₀=$v", dpi=300)
end

SAVE && savefig("april01-trajB.png")

let sol = ste01
    Random.seed!(4); trajplot(sol, 0:0.01:50, 1; limit=5, lab="", ms=5, msw=2, l=2)
    plot!(xguide="", yguide="", margin=0*Plots.mm, aspect_ratio=1)
    plot!(xticks=[], yticks=[])
    annotate!(0, 0, Plots.text("  start", 11, :black, :left))
    v = round(sol.speed; sigdigits=2)
    i = round(sol.info; sigdigits=2)
    plot!(title = "i/Dᵣ=$i, v/v₀=$v", dpi=300)
end

SAVE && savefig("april01-trajC.png")

let sol = abs01
    Random.seed!(8); trajplot(sol, 0:0.01:50, 1; limit=5, lab="", ms=5, msw=2, c=:darkorange, l=2)
    plot!(xguide="", yguide="", margin=0*Plots.mm, aspect_ratio=1)
    plot!(xticks=[], yticks=[])
    annotate!(0, 0, Plots.text("  start", 11, :black, :left))
    v = round(sol.speed; sigdigits=2)
    i = round(sol.info; sigdigits=2)
    plot!(title = "i/Dᵣ=$i, v/v₀=$v", dpi=300)
end

SAVE && savefig("april01-trajD.png")


#===== Paper =====#

fig2E = let
    scatter([0], [0], m=:+, c=:black, lab="") # "start")

    # Random.seed!(1); trajplot!(rev01, 0:0.01:50, 3; limit=5, lab="reverse", m=:uptri, ms=5, msw=2)
    # Random.seed!(5); trajplot!(tum01, 0:0.01:50, 3; limit=5, lab="tumble", ms=5, msw=2)
    # Random.seed!(2); trajplot!(abs01, 0:0.01:50, 3; limit=5, lab="steer", m=:downtri, ms=5, msw=2)

    Random.seed!(4); trajplot!(ste01, 0:0.01:50, 1; limit=5, lab="steer, odd", m=:diamond, ms=5, msw=2)
    Random.seed!(9); trajplot!(rev01, 0:0.01:50, 1; limit=5, lab="reverse", m=:uptri, ms=5, msw=2)
    # Random.seed!(6); trajplot!(tum01, 0:0.01:50, 1; limit=4.8, lab="tumble", ms=5, msw=2)
    Random.seed!(9); trajplot!(tum01, 0:0.01:50, 1; limit=5, lab="tumble", ms=5, msw=2)
    Random.seed!(2); trajplot!(abs01, 0:0.01:50, 1; limit=5, lab="steer, even", m=:downtri, ms=5, msw=2)


    scatter!([0], [0], m=:+, ms=20, msw=2, c=:black, lab="")
    annotate!(1.2, -0.5, Plots.text("start", 10, :black))

    # plot!(xguide="x", yguide="y", title="i/Dᵣ≈0.1, 0≤t≤50/Dᵣ")
    plot!(title="i/Dᵣ=0.1", xticks=[-5, 0, 5])
    plot!(legend=:topright, legendfontsize=10)
    plot!(legend=false)
    # plot!(size=(300,600), margin=Plots.mm*0, dpi=300, )
    plot!(size=(270,500), margin=Plots.mm*0, dpi=300, )
    # plot!(size=(300,600), margin=Plots.mm*0, dpi=300, )
end

# # Not sure the best format here...
# # PDF does something ugly to heatmaps,
# # and I need to combine these to add ABC letters anyway

SAVE && savefig(fig2E, "april01-fig2E.svg")




####################
##### FIGURE 3 #####


# pair1 = max2Wpair(0.5, 36; log=true)
# pair2 = max2Wpair(0.05, 72; log=true)
# pair3 = max2Wpair(0.005, 72; log=true, iterations=3000, g_abstol=1e-10)

function pairplot2(sol, left=false, bot=false)
  # color=cgrad([:lightgrey, :darkcyan, :darkblue])
  color=cgrad([:lightgrey, :darkblue, :darkgreen])
  # color=cgrad([:lightgrey, :black, :darkgreen])
  # color=cgrad([:lightgrey, :darkblue, :black])
  p1 = whole2heatmap(sol; fun=sqrt, both=false, color)
  plot!(p1; size=(300, 300), colorbar_title="√λ(Δθ,θ)")
  # plot!(p1; title = )
  left ? plot!(p1, yguide="turn angle Δθ") : plot!(p1, yguide="Δθ")
  bot ? plot!(p1, xguide="heading θ") : plot!(p1, xguide="θ")
  plot!(right_margin=Plots.mm*-5)

  dx = range(0,2pi,sol.n+1)[2:end-1]
  p2 = plot(dx, contactfun(sol)[2], xguide="Δθ", #, yguide="q̄(Δθ)                      ",
    c=:darkcyan, l=2, fill=0, fillalpha=0.1, lab="", yguidefontcolor=:darkcyan)
  plot!(p2, dx, contactfun(sol)[1], l=1, c=:darkmagenta, lab="")
  plot!(p2, size=(250, 250), yticks=[0,1], xticks=([0,pi,2pi], ["0", "π", "2π"]), xlims=[0,2pi])
  # left ? plot!(p2, xguide="turn angle Δθ") : plot!(p1, xguide="Δθ")
  plot!(p2, ylims=[0, 1.2])
  annotate!(p2, 5.5-left, 0.35, Plots.text("q̄(Δθ)", 10, :darkcyan))
  annotate!(p2, 5.5, 1.1, Plots.text("Ψ(Δθ)", 10, :darkmagenta))
  # plot!(p2, left_margin=0*Plots.mm)
  plot!(p2, left_margin=2*Plots.mm)

  # plot(p1, p2, size=(550, 300), layout=grid(1, 2, widths=normalize!([1.4, 1],1)))
  plot(p1, p2, size=(580, 300), layout=grid(1, 2, widths=normalize!([1.3, 1],1)))
  # p0 = plot(legend=false,grid=false,framestyle = :none,margin=0*Plots.mm) # axes=false,foreground_color_subplot=:white)
  # plot(p1, p0, p2; size=(520, 300), layout=@layout [a{0.58w} [_{0.05h}; b]])
end

fig3A = pairplot2(pair1, true)
fig3B = pairplot2(pair2)
fig3C = pairplot2(pair3, false, true)

SAVE && savefig(fig3A, "april01-fig3A.svg")
SAVE && savefig(fig3B, "april01-fig3B.svg")
SAVE && savefig(fig3C, "april01-fig3C.svg")



#===== Talk =====#


fig3nopsi = let
    p1 = whole2heatmap(pair1; fun=sqrt, both=false, color=cgrad([:lightgrey, :darkcyan, :darkblue]))
    plot!(p1, xguide="heading θ", yguide="turn angle Δθ", right_margin=-5*Plots.mm, left_margin=4*Plots.mm)

    p2 = whole2heatmap(pair2; fun=sqrt, both=false, color=cgrad([:lightgrey, :darkcyan, :darkblue]))
    plot!(p2, xguide="θ", yguide="", right_margin=-5*Plots.mm)

    p3 = whole2heatmap(pair3; fun=sqrt, both=false, color=cgrad([:lightgrey, :darkcyan, :darkblue]))
    plot!(p3, xguide="θ", yguide="", colorbar_title="√λ(Δθ,θ)", colorbar_fontcolor=:darkcyan)

    plot(p1, p2, p3,
        layout=(1,3), size=(900,250),
        bottom_margin=5*Plots.mm,
        dpi=300,
    )
end

SAVE && savefig(fig3nopsi, "april01-fig3nopsi.png")





#################################
##### FIGURE 4 -- 3D stuff ######


tumble3D = [max3sintumble(g) for g in logrange(0.0001, 1, 43)]
getfield.(tumble3D, :info)

reverse3D = [max3sinreverse(g) for g in logrange(0.0001, 1, 43)]
getfield.(reverse3D, :info)
reverse3Dstrong = [max3strongreverse(g, 100) for g in logrange(0.000001, 0.001, 23)]
getfield.(reverse3Dstrong, :info)

steer3D = [max3sindrift(g) for g in logrange(0.0001, 1, 33)]
getfield.(steer3D, :info)
steer3Dnew = [max3sindriftNEW(g) for g in logrange(0.00005, 3, 53)]
getfield.(steer3Dnew, :info)

let
    p = plot()
    myscatter!(p, steer3D, lab="old, D_c^2D = 1")
    plot!(p, x -> sqrt(x/6), 1e-2, 1; lab="old: sqrt(i/6)", s=:dash, c=:green)
    myscatter!(p, steer3Dnew, lab="new, D_c^1D = 2", c=:pink)
    plot!(p, x -> sqrt(x/3.5), 1e-2, 1; lab="new: sqrt(i/3.5)", c=:black)
    plot!(p, xaxis=:log10, yaxis=:log10)
end

flick3D = [max3sinflick(g) for g in logrange(0.002, 1, 23)]
getfield.(flick3D, :info)
getfield.(flick3D, :speed)
flick3Dstrong = [max3strongflick(g, 100) for g in logrange(0.0002, 0.002, 13)]  # I THINK THIS FAILED
getfield.(flick3Dstrong, :info)

# Wigner unsigned steering thing
roll3D = [max3wigner(@show(γ)) for γ in logrange(0.0005, 0.11, 23)]
roll3Dlo = [max3wigner(@show(γ), 0.1) for γ in logrange(0.0005, 0.12, 13)]
# roll3Dlolo = [max3wigner(@show(γ), 0.01) for γ in logrange(0.0002, 0.02, 13)]  # little change
roll3Dhi = [max3wigner(@show(γ), 10.0) for γ in logrange(0.0005, 0.04, 13)]
getfield.(roll3Dhi, :speed)
getfield.(roll3D, :info)

# roll3Dmore = [max3wigner(@show(γ)) for γ in range(0.0025, 0.025, 21)]
roll3Dlin = [max3wigner(@show(γ)) for γ in range(0.0004, 0.12, 53)]
roll3Dlolin = [max3wigner(@show(γ), 0.1) for γ in range(0.0004, 0.13, 53)]
roll3Dhilin = [max3wigner(@show(γ), 10.0) for γ in range(0.0004, 0.05, 53)]

noroll3D = [max3legendre(@show(γ)) for γ in logrange(0.001, 1, 13)]  # vibe-coded version of steering

# Optimal turn angles
sol_2s = max3sinjump(0.1, 16; iterations=1000)
sol_3s = max3sinjump(0.01, 32; iterations=1000)
# sol_4s = max3sinjump(0.001, 16, 32; iterations=1000, order=10)  # takes 40s
sol_4s = max3sinjump(0.001, 32, 32; iterations=1000, order=10)  # takes 3 mins!

fig4B = let
    maketitle(sol) = string("i/Dᵣ=", round(sol.info; sigdigits=2),
        ", v/v₀=", round(sol.speed; sigdigits=2))
    # color = cgrad([:lightgrey, :darkcyan, :darkblue])
    color = cgrad([:lightgrey, :darkblue, :black])

    p2 = jump3plot(sol_2s; color, fun=sqrt)
    plot!(p2, xguide="inclination θ", yguide="turn angle Δ", title=maketitle(sol_2s))

    p3 = jump3plot(sol_3s; color, fun=sqrt)
    plot!(p3, xguide="θ", yguide="", title=maketitle(sol_3s))

    p4 = jump3plot(sol_4s; color, fun=sqrt)
    plot!(p4, xguide="θ", yguide="", title=maketitle(sol_4s))

    plot(p2, p3, p4;
        colorbar=false,
        layout=(1,3),
        # size=(700,250),
        size=(660,230),
        bottom_margin=5*Plots.mm,
    )
end

SAVE && savefig(fig4B, "april01-fig4B.svg")


fig4A = let
    plot(xguide="information rate i/Dᵣ", yguide="climbing speed v/v₀", size=(550, 500))
    rotation = 40

    myplot!(steer3Dnew, lab="", fill=1, fillalpha=0.1, c=:white, l=0, fillcolor=:grey)
    # myscatter!(steer3D, lab="steering, OLD", m=:diamond, ms=3)
    myscatter!(steer3Dnew, lab="steering, direction ψ=0", m=:diamond, ms=3)
    # plot!(x -> sqrt(x/3), 1e-4, 1; lab="sqrt(i/3)", s=:dash, c=:red)

    # reverse
    c = COLOURS2D[:zeta]
    plot!(x -> sqrt(x/6), 1e-4, 1; lab="", s=:dash, c)
    plot!(x -> 1/2, 7, 1000; lab="", s=:dash, c)
    annotate!(0.19, 0.2, Plots.text("sqrt(i/6Dᵣ)", 9, c; rotation))
    annotate!(130, 0.45, Plots.text("v/v₀ → 1/2", 9, c))

    myscatter!(reverse3D, lab="reverse, Δ=π", m=:uptri)
    myscatter!(reverse3Dstrong, lab="", m=:square, alpha=0.5, ms=3)

    # # tumble
    # c=COLOURS2D[:lambda]
    # plot!(x -> sqrt(x/4), 1e-4, 1; lab="", s=:dash, c)
    # annotate!(0.1, 0.18, Plots.text("sqrt(i/4Dᵣ)", 9, c; rotation))

    myscatter!(tumble3D, lab="tumble")  # , all Δ

    # flick
    c = COLOURS2D[:kappa]
    plot!(x -> sqrt(x/12), 1e-4, 1; lab="", s=:dash, c)
    annotate!(0.22, 0.155, Plots.text("sqrt(i/12Dᵣ)", 9, c; rotation))

    myscatter!(flick3D, lab="flick, Δ=π/2", m=:pentagon)
    # myscatter!(flick3Dstrong, lab="", m=:square, alpha=0.5, ms=3)

    # optimal
    myscatter!([sol_2s, sol_3s, sol_4s]; c=:white, m=:star, lab="", ms=9, msw=0)
    myscatter!([sol_2s, sol_3s, sol_4s]; c=:black, m=:star, lab="optimal λ(Δ,θ)", ms=6)

    # roll
    myscatter!(roll3D, lab="steering, with Dψ=Dr", m=:downtri, c=:darkorange)
    myscatter!(roll3Dlin, lab="", m=:downtri, c=:darkorange)

    plot!(x -> x/8, 1e-4, 2; lab="", c=:darkorange)
    annotate!(0.17, 0.03, Plots.text("i/8Dᵣ", 9, :darkorange, rotation=55))

    # axes
    plot!(xaxis=:log10, yaxis=:log10)
    plot!(legendfontsize=9, legend=:bottomright)
    plot!(ylims=[0.017, 1.05], xlims=[0.09, 200])
    plot!(minorgrid=true)
    # plot!(size=(550, 380))
end

SAVE && savefig(fig4A, "april01-fig4A.pdf")




nothing

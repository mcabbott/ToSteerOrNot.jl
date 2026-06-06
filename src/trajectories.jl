export trajplot, trajplot!, trajectory


#####
##### trajplot, etc
#####

"""

```

trajplot(max2Wdrift(1.0; signed=true), 0:0.01:100, 5, lab="signed")
trajplot!(max2Wdrift(1.0; signed=false), 0:0.01:100, 5, lab="unsigned")

trajplot(max2Wreverse(1.0; log=true), 0:0.01:100, 5, lab="rverse")

trajplot(max2Wtumble(1.0; log=true), 0:0.01:100, 5, lab="tumble")

plot!(size=(400,800), aspect_ratio=:equal)

```
"""
function trajplot(nt::NamedTuple, times::AbstractRange, lines::Int=1; size=(400,800), kw...)
    plot(; xguide="x", yguide="y", size)
    trajplot!(nt, times, lines; kw...)
    # plot!([0], [0], m=:+, c=:black, ms=20, lab="", aspect_ratio=:equal) # start")
end

#=
Colors.parse(Colorant, :darkgreen) == RGB{N0f8}(0.0, 0.392, 0.0)
Colors.parse(Colorant, :green) == RGB{N0f8}(0.0, 0.502, 0.0)
=#

function trajplot!(nt::NamedTuple, times::AbstractRange, lines::Int=1; x0=0, limit=10, alpha=1, lab="", m=:circle, kw...)
    # Work out colours
    if hasproperty(nt, :mu)
        if nt.signed
            c, darkc = :red, :darkred
        else
            c, darkc = :orange, :darkorange
        end
     elseif hasproperty(nt, :jump)
        c, darkc = :grey, :darkgrey
    elseif hasproperty(nt, :lambda)
        c, darkc = :blue, :darkblue
    elseif hasproperty(nt, :zeta)
        # c, darkc = :lightgreen, :darkgreen
        c, darkc = Plots.RGB(0,0.65,0), Plots.RGB(0.1,0.4,0.1)
    else
        # error("nope!")
        c, darkc = :grey, :darkgrey
    end
    for line in 1:lines

        # Perform the simulation!
        xys_raw = trajectory(nt, times)
        xys_raw = map(xys_raw) do (x,y,z)
            (x+x0, y, z)
        end

        # Wrap to -10<x<10, inserting (NaN,NaN,true), (NaN,NaN,false)
        xys = trajlimit(xys_raw, limit)

        if hasproperty(nt, :jump)
            plot!(first.(xys), second.(xys); lab="", c, kw...)

            _isrev(θ) = 3.1 < θ < 3.2
            _isflick(θ) = (1.5 < θ < 4.8) && !_isrev(θ)
            _isother(θ) = (θ>0) && !_isrev(θ) && !_isflick(θ)

            xys3 = filter(t -> _isrev(t[3]), xys)  # reversals
            scatter!(first.(xys3), second.(xys3); lab="", c=:green, msw=0, ms=5, m=:uptri)

            xys3 = filter(t -> _isflick(t[3]), xys)  # flicks
            scatter!(first.(xys3), second.(xys3); lab="", c=:purple, msw=0, ms=4, m=:pentagon)

            xys3 = filter(t -> _isother(t[3]), xys)  # other
            scatter!(first.(xys3), second.(xys3); lab="", c=:orange, msw=0, ms=3)

        elseif hasproperty(nt, :lambda)
            plot!(first.(xys), second.(xys); lab="", c, kw...)

            xys3 = filter(t -> t[3], xys)  # points where it turns
            scatter!(first.(xys3), second.(xys3); lab="", c=:darkblue, msw=0, ms=2)
        else
            # Plot with alternating colours
            xys1 = filter(t -> t[3], xys)  # odd
            plot!(first.(xys1), second.(xys1); lab="", c, alpha, kw...)

            xys2 = filter(t -> !t[3], xys)  # even
            plot!(first.(xys2), second.(xys2); lab="", c=darkc, alpha, kw...)
        end

        plot!(first.(xys[end:end]), second.(xys[end:end]); lab=line==1 ? lab : "", m, c, kw...)
    end
    # hline!([nt.speed * maximum(times)]; c, s=:dash, lab="")
    plot!([x0], [0], m=:+, c=:black, ms=10, lab="", aspect_ratio=:equal)
end

mod3(x::Real, lo::Real, hi::Real) = mod(x-lo, hi-lo) + lo
mod3(xy::Tuple, lo::Real, hi::Real) = (mod3(xy[1], lo, hi), Base.tail(xy)...)

"""
    trajlimit([(x,y,b), (x,y,b), ...], 10)

This wraps to `-10 < x < 10`, inserting ponts with `NaN` to break the plot line.
"""
function trajlimit(xys::Vector, limit::Real)
    tmp = mod3.(xys, -limit, limit)
    diffx = abs.(first.(tmp[2:end]) .- first.(tmp[1:end-1]))
    inds = findall(>(limit), diffx)
    k = 1
    for i in inds
        insert!(tmp, i+k, (NaN, NaN, true))
        insert!(tmp, i+k+1, (NaN, NaN, false))
        k += 2
    end
    tmp
end

function trajectory(nt::NamedTuple, times::AbstractRange)
    if hasproperty(nt, :mmax)  # new 2026, see memory.jl
        xys = mem2traj(nt, times)
    elseif hasproperty(nt, :fun_raw)  # new 2026, see memory.jl
        xys = mem1traj(nt, times)
    elseif hasproperty(nt, :jump)
        xys = jumptraj(nt, times)
    elseif hasproperty(nt, :mu)
        xys = steertraj(nt, times)
    elseif hasproperty(nt, :lambda)
        xys = tumtraj(nt, times)
    elseif hasproperty(nt, :zeta)
        xys = revtraj(nt, times)
    else
        error("nope!")
    end
end

#=

zz1 = tuple.(1:100.0, 1:100.0, vcat(trues(43), falses(57)))
zz2 = SwimInfo.trajlimit(zz1, 10)
zz3 = vcat(filter(t -> t[3], zz2), filter(t -> !t[3], zz2))
plot(first.(zz3), getindex.(zz3,2))

=#



function tumtraj(nt::NamedTuple, times::AbstractRange)
    θ = 0.0  # current heading
    n = length(nt.lambda)  # nt.lambda is a vector of the tumble rates
    θs = angles2W(n)  # this is range(0, 2pi, length=n+1)[1:end-1]
    dθ = step(θs)  # = 2pi/n
    xy = (0.0, 0.0)
    traj = [(xy..., false)]  # will extend this vector with subsequent points
    dt = step(times)
    maximum(nt.lambda) * dt > 1 && @warn "rates λ are too high for this time-step" maxλ = maximum(nt.lambda) dt
    for t in times
        θ = mod2pi(θ)
        i = mod(searchsortedfirst(θs, θ-dθ/2), 1:n)  # approx index of heading θ
        if nt.lambda[i] * dt > rand()
            θ = rand(θs)
            # push!(traj, (NaN, NaN, false), (xy..., true))  # I forgot why I had NaN there
            push!(traj, (xy..., true))  # true means a tumble
        end
        θ += sqrt(2dt) * randn()  # sqrt(2 D dt) really, for diffusion
        xy = xy .+ sincos(θ) .* dt  # update position
        push!(traj, (xy..., false))
    end
    traj
end

#=

# sol_tum = max2Wtumble(1.0; log=true)

sol_tum = (objective = 0.061747132939097255, speed = 0.12173275529186295, info = 0.059985622352765694, gamma = 1.0, prob = [0.034601617772326974, 0.0344959012079262, 0.034182053519305186, 0.033669907101603606, 0.032975595338529326, 0.03212114367074013, 0.031133849040735895, 0.030045421371936585, 0.028890880668990634, 0.027707236695449534, 0.026532022316885496, 0.025401790959521105, 0.02435071436636976, 0.02340941312266148, 0.02260411718997138, 0.02195619272312061, 0.021482000801906237, 0.0211929931402405, 0.021095915755884615, 0.02119299314024052, 0.021482000801906275, 0.021956192723120677, 0.022604117189971455, 0.023409413122661565, 0.024350714366369847, 0.025401790959521216, 0.026532022316885617, 0.027707236695449656, 0.028890880668990752, 0.030045421371936696, 0.031133849040735975, 0.0321211436707402, 0.03297559533852936, 0.03366990710160362, 0.03418205351930519, 0.0344959012079262], tau = nothing, eta = true, lambda = [0.6057621486477696, 0.6107609312055187, 0.6258116477539292, 0.6510607503374858, 0.6866959432704143, 0.7328603722501954, 0.7895358089227195, 0.8564009572407253, 0.9326779200795836, 1.0169912682268192, 1.1072656140671242, 1.2006943991195416, 1.2938019876967082, 1.3826096931026992, 1.4628968122289798, 1.5305291755606103, 1.5818140708499335, 1.6138347996235132, 1.62472281267713, 1.613834799623514, 1.5818140708499355, 1.5305291755606105, 1.462896812228979, 1.382609693102699, 1.2938019876967075, 1.2006943991195405, 1.107265614067124, 1.016991268226819, 0.9326779200795837, 0.8564009572407258, 0.7895358089227197, 0.7328603722501956, 0.6866959432704145, 0.6510607503374854, 0.6258116477539294, 0.6107609312055188], positive = (lambda = log,), length = 36, calls = 144, time = 2.3)

function angles2W(n::Int)
    θs = range(0, 2pi, n+1)[1:end-1]
    @assert length(θs) == n
    θs
end

tumtraj(sol_tum, 0:0.01:1000)

=#

function steertraj(nt::NamedTuple, times::AbstractRange)
    θ = 0.0
    n = length(nt.mu)
    θs = angles2W(n)
    dθ = step(θs)
    xy = (0.0, 0.0)
    sig = true
    traj = [(xy..., sig)]
    dt = step(times)
    for t in times
        θ = mod2pi(θ)
        i = mod(searchsortedfirst(θs, θ-dθ/2), 1:n)
        θ += nt.mu[i] * dt
        if (nt.mu[i]>0) != sig  # then we change from left to right, or v-v
            sig = nt.mu[i]>0
            push!(traj, (NaN, NaN, sig), (xy..., sig))  # start a new segment, will find NaN when plotting
        end
        θ += sqrt(2*2*dt) * randn()  # sqrt(2 (Dr + Dc) dt) really, for diffusion
        xy = xy .+ sincos(θ) .* dt
        push!(traj, (xy..., sig))
    end
    traj
end

#=

# sol_signed = max2Wdrift(1.0; signed=true)

sol_signed = (objective = 0.11803033256038777, speed = 0.22684676601652531, info = 0.10881643345613755, gamma = 1.0, prob = [0.04161189153235155, 0.04134049812445779, 0.04053839927874255, 0.03926257895263635, 0.03759295463814197, 0.035624099691674016, 0.03345659465413531, 0.031189131722421787, 0.028912196815595367, 0.026703774770310655, 0.02462715477516414, 0.022730620677337717, 0.021048626344888233, 0.01960398059883848, 0.018410577724743075, 0.017476275848369117, 0.016805613440088966, 0.01640213706875325, 0.01627017462237215, 0.016270467329175184, 0.01668320127370011, 0.01736147977230471, 0.01830215021965656, 0.019501034623823563, 0.020950608118445284, 0.022637287543267874, 0.02453855555379272, 0.026620229030584495, 0.028834269104795397, 0.031117596981735447, 0.033392391875420056, 0.0355682724178663, 0.037546579060779524, 0.03922668290889975, 0.04051387803500448, 0.041328034869726075], mu = [-3.9898639947466563e-17, -0.07530856389630597, -0.22681531727278953, -0.37244568700752423, -0.5090272811319215, -0.6334130587629261, -0.7424883665403772, -0.8331929737235632, -0.9025639391433626, -0.9478063271289168, -0.9663993531600177, -0.9562446349543201, -0.9158597400870001, -0.8446129047846068, -0.7429828636905789, -0.612811626801608, -0.4575006615424715, -0.2820884129363437, -0.09314788895962729, 0.29047896090099645, 0.4656866964329066, 0.6206806836930844, 0.750450026755069, 0.8516232503883073, 0.9223869876378569, 0.9622872143186335, 0.9719751711074965, 0.952947331141456, 0.9073112221958861, 0.8375927990991184, 0.7465892314174568, 0.6372637333184459, 0.5126756747084923, 0.3759383755873348, 0.23019756134382738, 0.07862468147789137], signed = true, mu_raw = [-0.0016212005282234183, -0.07692976442452935, -0.2284365178010129, -0.3740668875357476, -0.5106484816601449, -0.6350342592911494, -0.7441095670686005, -0.8348141742517865, -0.904185139671586, -0.9494275276571401, -0.968020553688241, -0.9578658354825434, -0.9174809406152235, -0.8462341053128302, -0.7446040642188022, -0.6144328273298313, -0.4591218620706949, -0.2837096134645671, -0.09476908948785068, 0.28885776037277305, 0.4640654959046832, 0.619059483164861, 0.7488288262268457, 0.850002049860084, 0.9207657871096335, 0.9606660137904102, 0.9703539705792732, 0.9513261306132327, 0.9056900216676628, 0.8359715985708951, 0.7449680308892335, 0.6356425327902225, 0.511054474180269, 0.3743171750591114, 0.228576360815604, 0.07700348094966798], Dc = 0.9999999997584661, positive = (mu_raw = nothing, Dc = +), length = 37, calls = 1248, time = 1.6)

steertraj(sol_signed, 0:0.01:1000)

=#


function revtraj(nt::NamedTuple, times::AbstractRange)
    θ = 0.0
    n = length(nt.zeta)
    θs = angles2W(n)
    dθ = step(θs)
    xy = (0.0, 0.0)
    rev = false
    traj = [(xy..., rev)]
    dt = step(times)
    maximum(nt.zeta) * dt > 1 && @warn "rates are too high for this time-step" maxζ = maximum(nt.zeta) dt
    for t in times
        θ = mod2pi(θ)
        i = mod(searchsortedfirst(θs, θ-dθ/2), 1:n)
        if nt.zeta[i] * dt > rand()
            θ = mod2pi(pi + θ)
            rev = !rev
            push!(traj, (NaN, NaN, rev), (xy..., rev))  # start a new segment
        end
        θ += sqrt(2dt) * randn()  # sqrt(2 D dt) really, for diffusion
        xy = xy .+ sincos(θ) .* dt
        push!(traj, (xy..., rev))
    end
    traj
end


# mostly AI-written version for max2Wjump etc
function jumptraj(nt::NamedTuple, times::AbstractRange)
    Dr = 1.0
    v0 = 1.0

    jump = nt.jump
    nminus1, n = size(jump)
    nminus1 == n - 1 || error("jump should have size (n-1, n)")

    θs = angles2W(n)
    dθ = step(θs)

    θ = 0.0
    xy = (0.0, 0.0)

    traj = [(xy..., NaN)]

    dt = step(times)

    lambda = vec(sum(jump; dims=1))
    maximum(lambda) * dt > 1 && @warn "jump rates are too high for this time-step" maxλ = maximum(lambda) dt

    for t in times
        θ = mod2pi(θ)
        # approximate index of current heading θ
        i = mod(searchsortedfirst(θs, θ - dθ/2), 1:n)

        # Jump event
        if lambda[i] * dt > rand()

            # choose row r with probability jump[r,i] / λ[i]
            u = rand() * lambda[i]
            acc = 0.0
            rchosen = nminus1

            for r in 1:nminus1
                acc += jump[r, i]
                if u <= acc
                    rchosen = r
                    break
                end
            end
            turn = rchosen * dθ
            θ += turn
            push!(traj, (xy..., turn))
        end

        # Rotational diffusion
        θ += sqrt(2 * Dr * dt) * randn()
        # Swimming displacement
        xy = xy .+ v0 .* sincos(θ) .* dt

        push!(traj, (xy..., NaN))
    end
    traj
end

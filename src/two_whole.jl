export angles2W, whole2plot, jump2Wfold, whole2heatmap, signedinfo2Wrate, signedinfo2Wdrift, whole2diff, whole2box
export max2Wtumble, max2Wstrong, max2Wflick, max2Wstrongflick, max2Wdrift, max2Wreverse, max2Wsquare, max2Wleftright, max2Wstrongreverse
export max2Wjump, max2Wouter, max2Wmulti, max2Wfix, max2Wdiff, max2Wbox, max2Wreverseandflick, max2Wreverseorflick, max2Wreverseortumble, max2Wreversethenflick, max2Wupdownsteer, max2Wupdowntumble, max2Wupdownreverse
export steertojump, jumptosteer

#=

This file is based on "two_reflect.jl", which is my nice new cleanish thing...
but uses angles 0 to 2pi, stores the whole circle.

=#

function angles2W(n::Int)
    θs = range(0, 2pi, n+1)[1:end-1]
    @assert length(θs) == n
    θs
end

#####
##### Master equation matrices
#####

# IDK why these three are grouped... but they are!
function master2Wtrio(D::Real, lambda::Union{Nothing,Vector}, kappa::Union{Nothing,Vector}, mu::Union{Nothing,Vector}; delta::Int=0)
    n = length(something(lambda, kappa, mu))  # demands 1 non-nothing
    T = someeltype(lambda, kappa, mu)
    mat = zeros(T, n, n)

    master2Wdiff!(mat, D)
    master2Wtumble!(mat, lambda)
    master2Wflick!(mat, kappa, delta)
    master2Wdrift!(mat, mu)

    mat
end

function master2Wjump(D::Real, jump::Matrix) # , signed::Bool)
    n = size(jump, 1) + 1
    mat = zeros(eltype(jump), n, n)
    master2Wdiff!(mat, D)
    master2Wjump!(mat, jump) # , signed)
    mat
end

function master2Wreverse(D::Real, zeta::Vector)
    n = length(zeta)
    mat = zeros(eltype(zeta), n, n)
    master2Wdiff!(mat, D)
    master2Wreverse!(mat, zeta)
    mat
end

function master2Wdiff!(mat::Matrix, Diff::Real)
    n = size(mat,1)

    angles = angles2W(n)
    dθsquared = step(angles)^2

    for i in 1:n
        mat[i, i] -= 2*Diff/dθsquared
        mat[i, mod(i-1, 1:n)] += Diff/dθsquared
        mat[i, mod(i+1, 1:n)] += Diff/dθsquared
    end

    mat
end

##### Tumble

master2Wtumble!(mat::Matrix, lambda::Nothing, _...) = nothing
function master2Wtumble!(mat::Matrix, lambda::Vector)
    n = size(mat,1)
    for i in 1:n
        # sink:
        mat[i,i] -= lambda[i]
        # source -- uniform on the circle:
        for k in 1:n
            mat[k,i] += lambda[i] / n
        end
    end
    mat
end


##### Flick

master2Wflick!(mat::Matrix, kappa::Nothing, delta=0) = nothing
function master2Wflick!(mat::Matrix, kappa::Vector, delta::Int)
    delta > 0 || error("got delta = $delta, that's no good!")
    n = size(mat,1)
    for i in 1:n
        # sink:
        mat[i,i] -= kappa[i]
        # source -- equal prob left or right by 90º
        j1 = mod(i + delta, 1:n)
        j2 = mod(i - delta, 1:n)
        mat[i,j1] += kappa[j1] / 2
        mat[i,j2] += kappa[j2] / 2
    end
    mat
end

##### Reverse

function master2Wreverse!(mat::Matrix, rate::Vector)
    n = size(mat,1)
    for i in 1:n
        mat[i,i] -= rate[i]
        j = mod(i + n÷2, 1:n)
        mat[i,j] += rate[j]
    end
    mat
end

##### Drift == steering

master2Wdrift!(mat::Matrix, mu::Nothing) = nothing
function master2Wdrift!(mat::Matrix, mu::Vector)
    n = size(mat,1)
    @assert n == length(mu)
    dθ = step(angles2W(n))
    for i in 1:n
        drift = mu[i]/dθ
        if drift > 0
            ip1 = mod(i + 1, 1:n)
            mat[ip1,i] += drift
            mat[i,i] -= drift  # negative on the diagonal
        else
            im1 = mod(i - 1, 1:n)
            mat[i,i] += drift  # negative on the diagonal
            mat[im1,i] -= drift
        end
    end
    mat
end

##### Jump to anywhere

function master2Wjump!(mat::Matrix, jumprate::Matrix) #, signed::Bool)
    n = size(mat,1)
    # if signed
        @assert size(jumprate) == (n-1, n)  # jumprate[Δθ,θ]
        for i in 1:n  # θ
            for j in 1:n-1  # Δθ
                k = mod(i + j, 1:n)  # note that j=1 means shift by 1, there is no zero row
                mat[i, i] -= jumprate[j, i]  # sink
                mat[k, i] += jumprate[j, i]  # source
            end
        end
    # else
    #     @assert size(jumprate) == (n-1, n÷2+1)
    #     for i in axes(jumprate,2)  # θ
    #         i2 = n + 1 - i
    #         for j in 1:n-1  # Δθ
    #             k = mod(i + j, 1:n)  # note that j=1 means shift by 1, there is no zero row
    #             mat[i, i] -= jumprate[j, i]  # sink
    #             mat[i2, i2] -= jumprate[j, i]
    #             mat[k, i] += jumprate[j, i]  # source
    #             mat[k, i2] += jumprate[j, i]
    #         end
    #     end
    # end
    mat
end

function jump2Wfold(jumprate::Matrix)  # takes jumprate[Δθ,θ]
    n = size(jumprate, 1) + 1
    out = zeros(n, n)  # jump[θhat,θ]
    if size(jumprate) == (n-1, n)
        for i in 1:n
            for j in 1:n-1
                k = mod(i + j, 1:n)
                # out[i, i] -= jumprate[j, i]  # sink
                out[k, i] += jumprate[j, i]  # source
            end
        end
    # elseif size(jumprate) == (n-1, n÷2+1)
    #     for i in axes(jumprate,2)  # θ
    #         i2 = n + 1 - i
    #         for j in 1:n-1  # Δθ
    #             k = mod(i + j, 1:n)  # note that j=1 means shift by 1, there is no zero row
    #             # out[i, i] -= jumprate[j, i]  # sink
    #             # out[i2, i2] -= jumprate[j, i]
    #             out[k, i] += jumprate[j, i]  # source
    #             out[k, i2] += jumprate[j, i]
    #         end
    #     end
    end
    out
end

function fixjump2W(rate::Vector, dest::Vector)
    n = length(rate)
    @assert n == length(dest)
    out = zeros(eltype(rate), n-1, n)  # jumprate[Δθ,θ]
    for i in 1:n  # lamda[θ]
        for j in 1:n  # dest[θ+Δθ]
            k = mod(j-i, 1:n)
            k == n && continue
            out[k,i] += rate[i] * dest[j]
        end
    end
    out
end

#####
##### Speed & information
#####

function cos2W(prob::Vector)
    θs = angles2W(length(prob))
    sum(prob .* cos.(θs))
end

function sin2W(prob::Vector)
    θs = angles2W(length(prob))
    sum(prob .* sin.(θs))
end

"""
    eta2(p, λ, τ)

This is the prefactor ``1/(1+τ<λ>)`` which multiplies climb rate
when we make tumbles take finite time `τ`.

For `τ::Nothing` it returns `true`.
"""
eta2(prob::Vector, lambda::Vector, tau::Real) = 1 / (1 + tau * dot(prob, lambda))
eta2(prob::Vector, lambda::Vector, tau::Nothing) = true


function eta2Wquadratic(prob::Vector, jump::Matrix, tau::Real)
    n = length(prob)
    θs = angles2W(n)
    Δθs = vcat(θs[2:(n÷2+1)], θs[(n÷2+2):end] .- 2pi)
    1 / (1 + dot(tau .* Δθs.^2, jump, prob))
end
eta2Wquadratic(prob::Vector, jump::Matrix, tau::Nothing) = true

## This function is in "two_reflect.jl" but we want the same thing here:
# info2rate(prob::Vector, rate::Nothing) = false
# function info2rate(prob::AbstractVector, rate::Vector)
#     expect = dot(prob, rate)
#     sum(prob .* plogpoverq.(rate, expect))
# end

## And in "twodims.jl"
# plogpoverq(p::Real, q::Real) = (p*q<=0) ? zero(p*q) : p * log(p/q)

# info2drift(prob::Vector, Diff::Real, mu::Nothing) = false
# function info2drift(prob::Vector, Diff::Real, mu::Vector)
#     mumu = mean(mu)
#     varmu = sum(prob .* (mu .- mumu).^2)
#     varmu/Diff
# end

# function info2jump(prob::Vector, jumprate::Matrix)
#     pdelta = jumprate * prob
#     # jumprate * prob * log(jumprate / pdelta)
#     sum(@. prob' * plogpoverq(jumprate, pdelta))
# end

# function info2Wjump(prob::Vector, jumprate::Matrix, signed::Bool)
#     if signed
#         info2jump(prob, jumprate)
#     else
#         n = length(prob)
#         # @show size(prob) size(jumprate)
#         jump2 = hcat(jumprate, jumprate[:, reverse(2:end-1)])
#         # @show size(jump2)
#         info2jump(prob, jump2)
#     end
# end


function _full2Wjump(jumprate)
    n = size(jumprate, 1) + 1
    if size(jumprate) == (n-1, n)
        jumprate
    else
        hcat(jumprate, jumprate[:, reverse(2:end-1)])
    end
end

"""
    signedinfo2Wrate(prob, rate)

This is some attempt to separate `info2rate(prob, rate)` into `I(dθ, abs(θ)), I(dθ, sign(θ))`.
But it's a bit stupid, because the tuble / flick etc cases with rate always have symmetric prob anyway.

The drift case has asymmetric solutions.
"""
function signedinfo2Wrate(prob::Vector, rate::Vector)
    n = length(prob)
    iseven(n) || error("this doesn't handle odd n")

    p_sym = (prob .+ reverse(prob))./2
    Iabs = info2rate(p_sym, rate)

    angles = range(0, 2pi, n+1)[1:n]
    w1 = ((sin.(angles) .> eps()) .+ (sin.(angles) .> -eps())) ./ 2
    w2 = ((sin.(angles) .< eps()) .+ (sin.(angles) .< -eps())) ./ 2
    @assert all(≈(1), w1 .+ w2)

    p_sign = [dot(prob, w1), dot(prob, w2)]
    r_sign = [dot(rate, w1), dot(rate, w2)]
    @show p_sign r_sign
    Isign = info2rate(p_sign, r_sign)

    (; Iabs, Isign)
end


"""
    signedinfo2Wdrift(prob, rate)

This is some attempt to separate `info2drift(prob, D, mu)` into `I(dθ, abs(θ)), I(dθ, sign(θ))`.

But probably not correct!
"""
function signedinfo2Wdrift(prob::Vector, D::Real, mu::Vector)
    n = length(prob)
    iseven(n) || error("this doesn't handle odd n")

    p_sym = (prob .+ reverse(prob))./2
    Iabs = info2drift(p_sym, D, mu)  # or should it be abs.(mu)?

    angles = range(0, 2pi, n+1)[1:n]
    w1 = ((sin.(angles) .> eps()) .+ (sin.(angles) .> -eps())) ./ 2
    w2 = ((sin.(angles) .< eps()) .+ (sin.(angles) .< -eps())) ./ 2
    @assert all(≈(1), w1 .+ w2)

    p_sign = [dot(prob, w1), dot(prob, w2)]
    mu_sign = [dot(mu, w1), dot(mu, w2)]
    Isign = info2drift(p_sign, D, mu_sign)

    (; Iabs, Isign)
end


# function info2drift(prob::Vector, Diff::Real, mu::Vector, reflect::Bool=false)
#     if reflect  # the other half has opposite sign, hence mean is zero
#         varmu = sum(prob .* mu.^2)
#     else  # we have the whole circle, or something equivalent
#         mumu = mean(mu)
#         varmu = sum(prob .* (mu .- mumu).^2)
#     end
#     varmu/Diff
# end




#####
##### Optimisation -- simple strategies
#####

"""
    max2Wtumble(γ, n=36)

New version a lot like `optim2d(γ)`.
Fixes ``D_r=1``.
Keyword `tau` is time penalty per tumble, in units of ``1/D_r``.

```
sol = max2Wtumble(0.1; log=true, tau=0.5)
whole2plot(sol)
```
"""
function max2Wtumble(gamma::Real, n::Int=36; init=ones, tau=nothing, kw...)
    lambda = init(n)./10
    maximise((; lambda); kw...) do nt
        mat = master2Wtrio(1, nt.lambda, nothing, nothing)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)  # not sure this is faster here?
        eta = eta2(prob, nt.lambda, tau)
        speed = cos2W(prob) * eta
        info = info2rate(prob, nt.lambda)
        (; objective = speed - gamma * info, speed, info, gamma, prob, tau, eta)
    end
end

"""
    max2Wstrong(γ, n=100)

Tries to fit something like the Strong & Bialek solution, for tumbles alone.
Keyword `tau` is time penalty per tumble, in units of ``1/D_r``.
```
sol = max2Wstrong(0.01)
whole2plot(sol)
```
"""
function max2Wstrong(gamma::Real, n::Int=100; tau=nothing, kw...)
    θcrit = 0.3
    λflat = 100.0
    θlist = angles2W(n)
    maximise((; θcrit, λflat); kw...) do nt
        θc = clamp(nt.θcrit, 0+0.1, pi-0.1)
        # lambda is 0 or large:
        i = clamp(searchsortedfirst(θlist, θc), 1, n÷2)
        lambda = fill(nt.λflat, n)
        lambda[1:i] .= 0
        lambda[end+1-i:end] .= 0
        # now do linear interpolation?
        λstep = nt.λflat * (θlist[i] - θc)/step(θlist)
        lambda[i+1] = λstep
        lambda[end-i] = λstep

        mat = master2Wtrio(1, lambda, nothing, nothing)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        eta = eta2(prob, lambda, tau)
        speed = cos2W(prob) * eta
        info = info2rate(prob, lambda)
        (; objective = speed - gamma * info, speed, info, gamma, prob, lambda, tau, eta)
    end
end


"""
    max2Wflick(γ, n=36)

New version a lot like `optim1flick(γ)`.
Fixes ``D_r=1``.

```
sol = max2Wflick(0.1; log=true)
whole2plot(sol)
```
"""
function max2Wflick(gamma::Real, n::Int=36; init=ones, delta=div(n, 4), kw...)
    kappa = init(n)./10
    maximise((; kappa); kw...) do nt
        mat = master2Wtrio(1, nothing, nt.kappa, nothing; delta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, nt.kappa)
        (; objective = speed - gamma * info, speed, info, gamma, prob)
    end
end

function max2Wstrongflick(gamma::Real, n::Int=100; delta=div(n, 4), kw...)
    θcrit = 0.3
    λflat = 100.0
    θlist = angles2W(n)
    maximise((; θcrit, λflat); kw...) do nt
        θc = clamp(nt.θcrit, 0+0.1, pi-0.1)
        # lambda is 0 or large:
        i = clamp(searchsortedfirst(θlist, θc), 1, n÷2)
        lambda = fill(nt.λflat, n)
        lambda[1:i] .= 0
        lambda[end+1-i:end] .= 0
        # now do linear interpolation?
        λstep = nt.λflat * (θlist[i] - θc)/step(θlist)
        lambda[i+1] = λstep
        lambda[end-i] = λstep

        kappa = lambda
        mat = master2Wtrio(1, nothing, kappa, nothing; delta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, kappa)
        (; objective = speed - gamma * info, speed, info, gamma, prob, kappa)
    end
end

# function max2Wtumbleflick(gamma::Real, n::Int=36; init=ones, delta=div(n, 4))
#     lambda = init(n)./10
#     kappa = init(n)./10
#     maximise((; lambda, kappa)) do nt
#         prob = zero_eigvec(master2Wtrio(1, nt.lambda, nt.kappa, nothing; delta))
#         speed = cos2W(prob)
#         info = info2rate(prob, nt.lambda) + info2rate(prob, nt.kappa)
#         (; objective = speed - gamma * info, speed, info, gamma, prob)
#     end
# end

"""
    max2Wdrift(γ, n=36; signed=false)

New version of `optim1mu(γ)`, for what we're now calling "steering".
Default with `signed=false` is like `sym=true`, allows `μ(|θ|)` only.
Fixes ``D_r=1``.
"""
function max2Wdrift(gamma::Real, n::Int=36; init=zeros, signed::Bool=false)
    if signed
        mu_raw = init(n)./10
    else
        iseven(n) || error("haven't thought about odd n here")
        # include north and south, but east-west will be mirror images
        mu_raw = init(n÷2 + 1)./10
        mu_raw[1:3] .= -0.01
        mu_raw[end-3:end] .= +0.02
        # mu_raw[end÷2] = 0.1
        # mu_raw[end] = -0.1  # break the symmetry?
    end
    maximise((; mu_raw, Dc=1.0)) do nt
        if signed
            mu = nt.mu_raw .- mean(nt.mu_raw)
        else
            mu = vcat(nt.mu_raw, reverse(@view nt.mu_raw[2:end-1]))
            @assert length(mu) == n
        end
        mat = master2Wtrio(1+nt.Dc, nothing, nothing, mu)  # Diffusion term D_r + D_c
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        speed = cos2W(prob)
        speedx = sin2W(prob)
        info = info2drift(prob, nt.Dc, mu)  # Info term D_c only
        (; objective = speed - gamma * info, speed, speedx, info, gamma, prob, mu, signed)
    end
end

_sinlike(n::Int) = sin.(range(0,2pi,n))

"""
    max2Wsquare(γ, n=100; param=3)

This is a bit lie `max2Wstrong(γ)` for `max2Wdrift(γ; signed=false)`...
that is, it fits a square-wave profile for the turn rate.

With `param=3` it's my version, independent μplus & μminus.
With `param=2` it's Matt's version, with them locked together.
He invented this for trying to study this case at high information, analytically.

Slightly lower performance than full soln.

```
sol = max2Wsquare(0.001, 20)
sol.mu
whole2plot(sol)
```
"""
function max2Wsquare(gamma::Real, n::Int=100; param::Int=3, kw...)
    param in (2,3) || error("bad param keyword, must be 2 or 3")
    θcrit = 0.3
    μplus = 10.0
    μminus = 100.0  # minus this, near θ=0
    θlist = angles2W(n)
    maximise((; θcrit, μplus, μminus, Dc=1.0); kw...) do nt
        θc = clamp(nt.θcrit, 0+0.1, pi-0.1)
        i = clamp(searchsortedfirst(θlist, θc), 1, n÷2)
        mu = fill(nt.μplus, n)
        μminus = param==3 ? nt.μminus : nt.μplus
        mu[1:i] .= -μminus
        mu[end+1-i:end] .= -μminus
        # now do linear interpolation?
        λstep = (nt.μplus + μminus) * (θlist[i] - θc)/step(θlist) - μminus
        mu[i+1] = λstep
        mu[end-i] = λstep

        mat = master2Wtrio(1+nt.Dc, nothing, nothing, mu)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, gamma, prob, mu, signed=false)
    end
end

# function max2Wsquare(gamma::Real, n::Int=100; kw...)
#     θcrit = 0.3
#     μplus = 1.0
#     μminus = -10.0  # near θ=0...
#     θlist = angles2W(n)
#     maximise((; θcrit, μplus, μminus, Dc=1.0); kw...) do nt
#         θc = clamp(nt.θcrit, 0+0.1, pi-0.1)
#         i = clamp(searchsortedfirst(θlist, θc), 1, n÷2)
#         mu = fill(nt.μplus, n)
#         mu[1:i] .= nt.μminus
#         mu[end+1-i:end] .= nt.μminus
#         # now do linear interpolation?
#         λstep = (nt.μplus - nt.μminus) * (θlist[i] - θc)/step(θlist) + nt.μminus
#         mu[i+1] = λstep
#         mu[end-i] = λstep

#         prob = zero_eigvec(master2Wtrio(1+nt.Dc, nothing, nothing, mu))
#         speed = cos2W(prob)
#         info = info2drift(prob, nt.Dc, mu)
#         (; objective = speed - gamma * info, speed, info, gamma, prob, mu, signed=false)
#     end
# end


"""
    max2Wreverse(γ, n=36)

This is like `max2Wflick(γ)`, but 180º not 90º.
"""
function max2Wreverse(gamma::Real, n::Int=36; init=ones, kw...)
    zeta = init(n)./10
    maximise((; zeta); kw...) do nt
        mat = master2Wreverse(1, nt.zeta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, nt.zeta)
        (; objective = speed - gamma * info, speed, info, gamma, prob)
    end
end

function max2Wstrongreverse(gamma::Real, n::Int=100; delta=div(n, 4), kw...)
    θcrit = 0.3
    λflat = 100.0
    θlist = angles2W(n)
    maximise((; θcrit, λflat); kw...) do nt
        θc = clamp(nt.θcrit, 0+0.1, pi-0.1)
        # lambda is 0 or large:
        i = clamp(searchsortedfirst(θlist, θc), 1, n÷2)
        lambda = fill(nt.λflat, n)
        lambda[1:i] .= 0
        lambda[end+1-i:end] .= 0
        # now do linear interpolation?
        λstep = nt.λflat * (θlist[i] - θc)/step(θlist)
        lambda[i+1] = λstep
        lambda[end-i] = λstep

        zeta = lambda
        mat = master2Wreverse(1, zeta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, zeta)
        (; objective = speed - gamma * info, speed, info, gamma, prob, zeta)
    end
end

#####
##### Optimisation -- mixed strategies
#####

"""
    max2Wreverseandflick(γ, n=36)

This uses the same rate for both reverse *and* flick.
That might simulate a bug which alternates reverse-flick-reverse, in steady state?
Not sure, this is one strategy for such bugs, maybe you could do better.
"""
function max2Wreverseandflick(gamma::Real, n::Int=36; init=ones, delta=div(n, 2), kw...)
    kappa = init(n)./10
    maximise((; kappa); kw...) do nt
        mat = master2Wtrio(1, nothing, nt.kappa, nothing; delta)
        master2Wreverse!(mat, nt.kappa)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, 2 .* nt.kappa)  # pay just one information... maybe you want a 2 though?
        (; objective = speed - gamma * info, speed, info, gamma, prob, zeta=nt.kappa)
    end
end

"""
    max2Wreverseorflick(γ, n=36)

This uses independent rates for reverse & flick.
That should simulate a bug which can choose which to perform?

NOT SURE THIS WORKS!
"""
function max2Wreverseorflick(gamma::Real, n::Int=36; init=ones, delta=div(n, 4), kw...)
    kappa = init(n)./10
    zeta = init(n)./10
    maximise((; kappa, zeta); kw...) do nt
        mat = master2Wtrio(1, nothing, nt.kappa, nothing; delta)
        master2Wreverse!(mat, nt.zeta)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, nt.kappa) + info2rate(prob, nt.zeta)  # two indep choices
        (; objective = speed - gamma * info, speed, info, gamma, prob)
    end
end


function max2Wreverseortumble(gamma::Real, n::Int=36; init=ones, kw...)
    lambda = init(n)./10
    zeta = init(n)./10
    maximise((; lambda, zeta); kw...) do nt
        mat = master2Wtrio(1, nt.lambda, nothing, nothing)
        master2Wreverse!(mat, nt.zeta)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, nt.lambda) + info2rate(prob, nt.zeta)  # two indep choices
        (; objective = speed - gamma * info, speed, info, gamma, prob)
    end
end


"""
    max2Wreversethenflick(γ, n=36)

This simulates reverse-then-flick by keeping independent probabilities for even and odd runs?

NOT SURE THIS IS CORRECTLY DONE

```
sol1 = max2Wreversethenflick(0.03; log=true)
whole2plot(sol1; lock=true)

sol3 = max2Wreversethenflick(0.3; log=true)
whole2plot(sol3; lock=true)
```
"""
function max2Wreversethenflick(gamma::Real, n::Int=36; init=ones, delta=div(n, 4), kw...)
    kappa = init(n)./10
    zeta = init(n)./10
    maximise((; kappa, zeta); kw...) do nt
        mat1 = master2Wtrio(1, nothing, nt.kappa, nothing; delta)  # flick
        mat2 = master2Wreverse(1, nt.zeta)  # reverse
        mat = [0I mat1; mat2 0I]  #
        prob12 = reshape(zero_eigvec_iter(mat),n,2)
        # prob12 = reshape(zero_eigvec_iter3(mat),n,2)
        prob1, prob2 = eachcol(prob12)
        prob = prob1 + prob2
        @assert sum(prob) ≈ 1
        speed = cos2W(prob)
        info = sum(prob1) .* info2rate(prob1 ./ sum(prob1), nt.zeta) + sum(prob2) .* info2rate(prob2 ./ sum(prob2), nt.kappa)   # MIGHT BE BACKWARDS
        (; objective = speed - gamma * info, speed, info, gamma, prob, prob1, prob2, mat)
    end
end

#####
##### Optimisation -- binary strategies
#####


"""
    max2Wleftright(γ, n=100)

This is like `max2Wdrift(γ; signed=true)` but fits a square wave, steering at constant `±μ` for left/right.
Simpler than `max2Wsquare(γ)` which is a square wave for `max2Wdrift(γ; signed=false)`,
as there's only one parameter, not three.
This is the case which Matt solved exactly, analtically.

```
sol = max2Wleftright(0.1, 20)
sol.mu
whole2plot(sol)
```
"""
function max2Wleftright(gamma::Real, n::Int=100; kw...)
    μplus = 0.1
    iseven(n) || error("haven't thought about odd n here")
    θlist = angles2W(n)
    @assert θlist[1] ≈ 0
    maximise((; μplus, Dc=1.0); kw...) do nt
        mu = fill(-nt.μplus, n)
        mu[1] = 0  # exactly uphill
        mu[n÷2+2:end] .= nt.μplus
        mu[n÷2+1] = 0
        mat = master2Wtrio(1+nt.Dc, nothing, nothing, mu)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, gamma, prob, mu, signed=true)
    end
end

"""
    max2Wupdowntumble(γ, n=100)

2-parameter ansatz for tumbles, which uses only whether it's going uphill or down.
(Derek's suggestion.)

```
sol = max2Wupdowntumble(1.0)
whole2plot(sol)
```
"""
function max2Wupdowntumble(gamma::Real, n::Int=100; kw...)
    λplus = 0.1
    λminus = 0.1
    mod(n,4)==0 || error("please make sure n is divisible by 4!")
    θlist = angles2W(n)
    @assert θlist[1] ≈ 0
    maximise((; λplus, λminus, Dc=1.0); kw...) do nt
        lambda = fill(nt.λplus, n)
        lambda[n÷4+2:3n÷4] .= nt.λminus
        lambda[n÷4+1] = (nt.λplus + nt.λminus)/2  # on the contour, 90 degrees
        lambda[3n÷4+1] = (nt.λplus + nt.λminus)/2  # ditto, 270 degrees
        mat = master2Wtrio(1, lambda, nothing, nothing)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, lambda)
        (; objective = speed - gamma * info, speed, info, gamma, prob, lambda)
    end
end

function max2Wupdownreverse(gamma::Real, n::Int=100; kw...)
    λplus = 0.1
    λminus = 0.1
    mod(n,4)==0 || error("please make sure n is divisible by 4!")
    θlist = angles2W(n)
    @assert θlist[1] ≈ 0
    maximise((; λplus, λminus, Dc=1.0); kw...) do nt
        zeta = fill(nt.λplus, n)
        zeta[n÷4+2:3n÷4] .= nt.λminus
        zeta[n÷4+1] = (nt.λplus + nt.λminus)/2  # on the contour, 90 degrees
        zeta[3n÷4+1] = (nt.λplus + nt.λminus)/2  # ditto, 270 degrees
        mat = master2Wreverse(1, zeta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2rate(prob, zeta)
        (; objective = speed - gamma * info, speed, info, gamma, prob, zeta)
    end
end


#####
##### Optimisation -- complex strategies
#####


"""
    max2Wjump(γ, n=20; signed=true, mirror=false)

This optimises a matrix `size(jump) == (n-1,n)` representing rates `λ(Δθ|θ)`
for any nonzero `Δθ`.

Case `signed=false` has matrix `size(jump) == (n-1,n÷2+1)`,
constraining the rates to depend only on `abs(θ)`.
But `Δθ` is free -- it can still decide to turn left or right.

Option `mirror=true` ensures that `Δθ` and `-Δθ` have the same weight.

Keyword `tau2` is the coefficient of a time delay ``τ₂ Δθ²``, default `nothing` meaning zero.
"""
function max2Wjump(gamma::Real, n::Int=20; signed::Bool=true, mirror::Bool=false, tau2=nothing, init=ones, kw...)
    if signed && mirror
        # error("signed=true, mirror=false isn't supported right now")
        iseven(n) || error("haven't thought about odd n here")
        jump_raw = init(n÷2, n)./10
    elseif signed
        @assert !mirror
        jump_raw = init(n-1, n)./10
    elseif mirror
        @assert !signed
        iseven(n) || error("haven't thought about odd n here")
        jump_raw = init(n÷2, n÷2 + 1)./10
    else
        @assert !signed
        @assert !mirror
        iseven(n) || error("haven't thought about odd n here")
        jump_raw = init(n-1, n÷2 + 1)./10
    end
    maximise(_max2Wjump(; gamma, n, signed, mirror, tau2), (; jump_raw); kw...)
end
function _max2Wjump(; gamma, n, signed, mirror, tau2, kw...)
    function objective(nt::NamedTuple)
        jump = if signed && mirror
            vcat(nt.jump_raw, @view nt.jump_raw[reverse(1:end-1), reverse(1:end)])
        elseif signed
            @assert !mirror
            nt.jump_raw
        elseif mirror
            @assert !signed
            tmp = vcat(nt.jump_raw, @view nt.jump_raw[reverse(1:end-1), :])
            hcat(tmp, @view tmp[:, reverse(2:end-1)])
        else
            @assert !mirror && !signed
            hcat(nt.jump_raw, @view nt.jump_raw[:, reverse(2:end-1)])
        end
        mat = master2Wjump(1, jump)
        # prob = zero_eigvec_iter(mat)
        # prob = zero_eigvec_iter3(master2Wjump(1, jump))
        prob = zero_eigvec_iter3(mat)
        eta = eta2Wquadratic(prob, jump, tau2)
        speed = cos2W(prob) * eta
        info = info2jump(prob, jump)
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, n, tau2, eta, signed, mirror)
    end
end

"""
    max2Wouter(γ, n=20; mirror=false)

This optimises a matrix `jump = que * lambda'` representing rates ``λ(Δθ|θ)``,
factorised into ``λ(θ)`` and ``q(Δθ)``.

With `mirror=false`, it finds ``q(Δθ) = δ(Δθ+α)`` with nonzero offset `α` at high info,
in which case ``λ(θ)`` also ends up asymmetric.

With `mirror=true`, it imposes that ``q(Δθ)`` be even.

```
sol1 = max2Wouter(0.1)
sol2 = max2Wouter(0.01)
sol3 = max2Wouter(0.01; mirror=true)
whole2heatmap([sol1, sol2, sol3])

info2jump(sol1.prob, sol1.jump)
info2rate(sol1.prob, sol1.lambda)
```
"""
function max2Wouter(gamma::Real, n::Int=20; mirror::Bool=false, tau2=nothing, init=ones, kw...)
    lambda = init(n)./10
    if mirror
        iseven(n) || error("haven't thought about odd n here")
        que_raw = -init(n÷2)./10
    else
        # Making this -ve is a hack to make log=true act on lambda alone
        que_raw = -init(n-1)./10
    end
    maximise(_max2Wouter(; gamma, n, mirror, tau2), (; lambda, que_raw); kw...)
end
# Not sure I love this setup, but the point is to be able to call objective later, outside of maximise.
function _max2Wouter(; gamma, n, mirror, tau2, kw...)
    function objective(nt::NamedTuple{(:lambda, :que_raw)})
        _raw = @. clamp(-nt.que_raw, 0, Inf)
        _que = if mirror
            vcat(_raw, @view _raw[reverse(1:end-1)])
        else
            _raw
        end
        que = _que ./ sum(_que)
        jump = que .* transpose(nt.lambda)
        @assert size(jump) == (n-1, n)
        mat = master2Wjump(1, jump)
        # prob = zero_eigvec_iter(mat)
        prob = zero_eigvec_iter3(mat)
        eta = eta2Wquadratic(prob, jump, tau2)
        speed = cos2W(prob) * eta
        # info = info2jump(prob, jump)  # what I first wrote. Finds solutions as in draft -04, with bifurcations
        info = info2rate(prob, nt.lambda)  # this... now also works!
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, que, n, tau2, eta, mirror)
    end
end

"""
    max2Wmulti(γ, n=20)

This works out jump rate as ``λ(Δθ|θ) = ∑ᵢ λᵢ(θ) qᵢ(Δθ)``, with 3 terms.
Each term is like `max2Wouter(; mirror=true)`
but crucially the information is like `max2Wjump` i.e. it has cross-terms between the `qᵢ(Δθ)`.

```
sols = [max2Wmulti(γ, 36; init=ones) for γ in [0.3, 0.1, 0.03]];
whole2heatmap(sols; size=(1000,400), layout=(1,3))

sols_r = [max2Wmulti(γ, 36; init=rand) for γ in [0.3, 0.1, 0.03]];
whole2heatmap(sols_r; size=(1000,400), layout=(1,3))
```
"""
function max2Wmulti(gamma::Real, n::Int=20; init=ones, kw...)
    iseven(n) || error("haven't thought about odd n here")
    lambda_raw = init(n÷2 + 1, 4)./10
    que_raw = init(n÷2, 4)./10

    # Try to prime it?
    que_raw[n÷4 : end,1] .= 1  # reverse
    que_raw[n÷4,2] = 1
    que_raw[n÷3,3] = 1

    # que_raw[n÷4 : 3n÷4,1] = 1  # reverse
    # que_raw[n÷4,2] = que_raw[3n÷4,2] = 10
    # que_raw[n÷3,3] = que_raw[2n÷3,3] = 10

    maximise(_max2Wmulti(; gamma, n), (; lambda_raw, que_raw); kw...)
end
function _max2Wmulti(; gamma, n, kw...)
    function objective(nt::NamedTuple{(:lambda_raw, :que_raw)})
        lambdas = eachcol(vcat(nt.lambda_raw, @view nt.lambda_raw[reverse(2:end-1),:]))
        _que = vcat(nt.que_raw, @view nt.que_raw[reverse(1:end-1),:])
        ques = eachcol(_que ./ sum(_que))
        jump = (ques[1] .* transpose(lambdas[1])
                .+ ques[2] .* transpose(lambdas[2])
                .+ ques[3] .* transpose(lambdas[3])
                 .+ ques[4] .* transpose(lambdas[4]))
        @assert size(jump) == (n-1, n)
        mat = master2Wjump(1, jump)
        # prob = zero_eigvec_iter(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        # info = info2rate(prob, lambdas[1]) + info2rate(prob, lambdas[2]) + info2rate(prob, lambdas[3])  # No this is WRONG, you forgot the cross terms
        info = info2jump(prob, jump)
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, ques, lambdas, n, mirror=true, signed=false)
    end
end


"""
    max2Wfix(γ, n=20; mirror=false)

This optimises another factorised form, ``λ(Δθ|θ) = λ(θ) q(θ+Δθ)``, i.e. fixed output distribution of angles.
Ben's claim is that this is the (or at least an) optimal solution for general ``λ(Δθ|θ)`` case, i.e. using sign of ``θ`` and not imposing symmetry on ``Δθ``.

```
sol1 = max2Wfix(0.1; log=true)  # this function
sol1.prob
sol1.qout  # output angle distribution

sol2 = max2Wjump(0.1)  # as before

plot(whole2heatmap(sol1), whole2heatmap(sol2), size=(800, 800))
```
"""
function max2Wfix(gamma::Real, n::Int=20; mirror::Bool=false, init=ones, kw...)
    lambda = init(n)./10
    que_raw = init(n)./10
    maximise((; lambda, que_raw); kw...) do nt
        qout = nt.que_raw ./ sum(nt.que_raw)
        jump = fixjump2W(nt.lambda, qout)
        @assert size(jump) == (n-1, n)
        mat = master2Wjump(1, jump)
        # prob = zero_eigvec_iter(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2jump(prob, jump)  # NB this full one!
        # info = info2rate(prob, nt.lambda)  # this is different, soln is a delta-function in q.
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, qout, n, mirror=false, signed=true)
    end
end

"""
    max2Wdiff(γ, n=20, m=20; Dt=10, control=true)

Here the bug gets to decide how long to spend turning by diffusion with ``Dt >> 1``.
It controls the rate ``λ(τ|θ)`` for starting tumbles,
which end deterministically time ``τ`` later.

Allows `m` times ``τ`` from `1/10Dt` to `5/Dt` (log spaced).
That's about the right range of tumble times?
Calling `heatmap2Wdiff(n, m; Dt)` will draw the matrix ``ρ(Δθ|τ)`` being used.

`control=true` uses information from θ to τ.
`control=false` uses information from θ to Δθ instead.

```
sol = max2Wdiff(0.01; Dt=100)
whole2diff(sol)

sol2 = max2Wdiff(0.01; Dt=100, control=false)
whole2diff(sol2)
```

Ben had some complaint...
"""
function max2Wdiff(gamma::Real, n::Int=20, m::Int=20; Dt::Real=10, control::Bool=true, init=ones, kw...)
    rates = init(m, n)./10  # this is λ(τ|θ)
    Δθs = angles2W(n)[2:end]
    @assert !(0 in Δθs)

    taus = logrange(1/(10*Dt), 5/Dt, length=m)
    _rho = [circdiff(Δθ, t, Dt) for Δθ in Δθs, t in taus]
    rho = _rho ./ sum(_rho; dims=1);  # normalise over θ
    # heatmap(rho)  # this result is independent of Dt
    @assert size(rho,1) == n-1
    @assert all(isfinite, rho)

    maximise((; rates); kw...) do nt
        jump = rho * nt.rates  # this is λ(Δθ|θ)
        # prob = zero_eigvec_iter(master2Wjump(1, jump))
        prob = zero_eigvec_iter3(master2Wjump(1, jump))
        expect = sum(taus .* nt.rates .* transpose(prob))  # ∫dτ τ ∫dθ p(θ) λ(τ|θ)
        eta =  1/ (1 + expect)
        speed = cos2W(prob) * eta
        info = info2jump(prob, control ? nt.rates : jump)
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, taus, Dt, control, eta)
    end
end

# 1/(2 Pi) + 1/Pi Sum[Exp[-Dr  n^2  t] Cos[n  x], {n, 0, 100}]
function circdiff(Δθ::Real, t::Real, D::Real; order=20)
    p = 1/2 + sum(exp(-D*k^2*t)*cos(k*Δθ) for k in 0:order)
    # t == 0 ? oftype(p, θ == 0) : p  # not really needed
end

export heatmap2Wdiff

function heatmap2Wdiff(n::Int=35, m::Int=20; Dt::Real=10)
    # THIS CODE ALL COPIED VERBATIM FROM JUST ABOVE!
    Δθs = angles2W(n)[2:end]
    @assert !(0 in Δθs)

    taus = logrange(1/(10*Dt), 5/Dt, length=m)
    _rho = [circdiff(Δθ, t, Dt) for Δθ in Δθs, t in taus]
    rho = _rho ./ sum(_rho; dims=1);  # normalise over θ

    # THIS IS NEW
    heatmap(taus, Δθs, rho;
        yguide="change Δθ (n=$n)", yticks=([0, π, 2pi], ["0", "π", "2π"]),
        xguide="delay τ (m=$m)", xaxis=:log10,
        title="matrix ρ(Δθ|τ) for Dt/Dr=$Dt",
    )
end


"""
    max2Wbox(γ, n=20)

Here the bug gets to decide what the maximum turn angle is.
That was exactly my march meeting setup.

```
sol = max2Wbox(0.01)
whole2box(sol)
```

TODO make this take `tau` which adds an incentive for smaller turns?
"""
function max2Wbox(gamma::Real, n::Int=20; init=ones, kw...)
    rates = init((n+1)÷2, n)./10
    _rho = [(j≤i || (n-j)≤i) for j in 1:(n-1), i in 1:((n+1)÷2)]
    rho = _rho ./ sum(_rho; dims=1)
    # @assert size(rho,1) == n-1
    # @show size(rho) size(rates)
    maximise(_max2Wbox(; gamma, n, rho), (; rates); kw...)
end
function _max2Wbox(; gamma, n, rho, kw...)
    @assert size(rho) == (n-1, (n+1)÷2)
    function objective(nt::NamedTuple)
        # @show size(rho) size(nt.rates)
        @assert size(nt.rates) == ((n+1)÷2, n)
        jump = rho * nt.rates
        mat = master2Wjump(1, jump)
        # prob = zero_eigvec(mat)
        # prob = zero_eigvec_iter(mat)
        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2jump(prob, nt.rates)  # control not jump?
        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, n, rho)
    end
end

#####
##### Conversion
#####

"""
    steer2jump(nt, δ=1)

Takes the NamedTuple returned by `max2Wdrift`, describing a steering solution,
and makes one like that from `max2Wjumnp`.

prob, speed, info are re-computed.
Prints out the comparison between old & new.

```
sol1 = max2Wdrift(0.3, 72; signed=true)
sol2 = steertojump(sol1, 2)
plot(whole2plot(sol1), whole2heatmap(sol2; both=false), plot_title="mu-to-jump")
```
"""
function steertojump(nt::NamedTuple, delta::Int=1)
    n = length(nt.prob)
    θs = angles2W(n)

    deltatheta = delta * step(θs)
    fplus = nt.Dc/deltatheta^2 .+ nt.mu ./ (2*deltatheta)
    fminus = nt.Dc/deltatheta^2 .- nt.mu ./ (2*deltatheta)
    if any(<(0), fplus) || any(<(0), fminus)
        @warn "generating negative jump rates, that's not good"
    end
    jump = zeros(n-1, n)
    jump[delta, :] .= clamp.(fplus, 0, Inf)
    jump[end+1-delta, :] .= clamp.(fminus, 0, Inf)

    prob = zero_eigvec_iter3(master2Wjump(1, jump))
    speed = cos2W(prob)
    info = info2jump(prob, jump)

    old = nt; @info "steer to jump" old.info info old.speed speed old.gamma

    (; speed, info, gamma=NaN, delta, deltatheta, fplus, fminus, jump, prob, mirror=false, signed=nt.signed)
end

function steertojump(nt::NamedTuple, deltas::AbstractVector{Int}, wei::AbstractVector=ones(length(deltas)))
    n = length(nt.prob)
    θs = angles2W(n)

    jump = zeros(n-1, n)

    scales = wei ./ sum(wei)

    for (scale, delta) in zip(scales, deltas)
        deltatheta = delta * step(θs)
        fplus = nt.Dc/deltatheta^2 .+ nt.mu ./ (2*deltatheta)
        fminus = nt.Dc/deltatheta^2 .- nt.mu ./ (2*deltatheta)
        if any(<(0), fplus) || any(<(0), fminus)
            @warn "generating negative jump rates, that's not good"
        end
        jump[delta, :] .= clamp.(fplus, 0, Inf) .* scale
        jump[end+1-delta, :] .= clamp.(fminus, 0, Inf) .* scale
    end

    prob = zero_eigvec_iter3(master2Wjump(1, jump))
    speed = cos2W(prob)
    info = info2jump(prob, jump)

    old = nt; @info "steer to jump" old.info info old.speed speed old.gamma

    (; speed, info, gamma=NaN, deltas, jump, prob, mirror=false, signed=nt.signed)
end


"""
```
sol3 = max2Wjump(0.3, 20, log=true)
sol4 = jumptosteer(sol3)

sol4.Dc
plot(whole2heatmap(sol3; both=false), whole2plot(sol4),  plot_title="jump-to-mu")
```
"""
function jumptosteer(nt::NamedTuple)
    n = length(nt.prob)
    θs = angles2W(n)
    Δθs = mod2pi.(θs[2:end] .+ pi) .- pi

    mu = zeros(n)
    Dc = 0.0
    for j in 1:n
        for k in 1:n-1
            mu[j] += Δθs[k] * nt.jump[k,j]
            Dc += Δθs[k]^2 * nt.jump[k,j] / (2n)
        end
    end

    mat = master2Wtrio(1+Dc, nothing, nothing, mu)  # Diffusion term D_r + D_c
    prob = zero_eigvec_iter3(mat; shift=1e-10)
    speed = cos2W(prob)
    info = info2drift(prob, Dc, mu)  # Info term D_c only

    old = nt; @info "jump to steer" old.info info old.speed speed old.gamma

    (; speed, info, gamma=NaN, mu, Dc, prob, signed=true)
end

#####
##### Plotting -- simple strategies, one function
#####

function whole2plot(nt::NamedTuple; yguide="", yguide2="", xguide="θ (n=$(length(nt.prob)))", lock::Bool=true, sigdigits=3, yticks=:auto)
    n = length(nt.prob)
    # theta = angles2W(n)
    theta = range(0, 2pi, n+1)  # wrap around
    prob = vcat(nt.prob, nt.prob[1])  # wrap around
    p = plot(theta, prob; c=:grey, fill=0, fillalpha=0.1,
        yguide, yticks=[], ylim=[0, 1.1*maximum(prob)],
        xguide, xlims=[0,2pi],
        xticks=([0, pi, 2pi], ["0","π","2π"]))
    if hasproperty(nt, :prob1)  #
        prob1 = vcat(nt.prob1, nt.prob1[1])
        plot!(p, theta, prob1; c=:grey, s=:dot)
        prob2 = vcat(nt.prob2, nt.prob2[1])
        plot!(p, theta, prob2; c=:grey, s=:dash)
    end
    cnt = 0
    for sy in keys(COLOURS2D)
        hasproperty(nt, sy) || continue
        c = COLOURS2D[sy]
        if sy === :mu && nt.signed === false
            c = COLOURS2D[:abs]  # dark orange for unsigned steering
        elseif sy === :mu && hasproperty(nt, :μplus)
            c = COLOURS2D[:sign]  # dark red for left-right steering
        end
        extra = lock ? _find_ylim(nt) : (;)
        if sy === :mu && nt.signed === true && nt.gamma > 1
            # Matt's analytic low-info signed steering
            plot!(twinx(), x -> -sin(x)*sqrt(8*nt.info), 0, 2pi; c, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi])
        elseif sy === :mu && nt.signed === false && nt.gamma > 0.1
            # Matt's analytic low-info unsigned steering
            plot!(twinx(), x -> 2-cos(x)*sqrt(8*nt.info), 0, 2pi; c, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi])
            # @show 2-cos(pi)*sqrt(8*nt.info)
        elseif sy === :lambda && nt.gamma > 1
            # Jose's solution for low-info tumbles:
            # plot!(twinx(), x -> 1 -    cos(x)*sqrt(8*nt.info), 0, 2pi; c, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi])
            plot!(twinx(), x -> 1 -    cos(x)*sqrt(4*nt.info), 0, 2pi; c, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi])  # better version, november
        elseif sy === :zeta && nt.gamma > 1
            # Just a guess?
            plot!(twinx(), x -> 0.5 - 0.5*cos(x)*sqrt(8*nt.info), 0, 2pi; c, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi], yticks)
            # plot!(twinx(), x -> 0.5 -  0.5*cos(x)*sqrt(4*nt.info), 0, 2pi; c, l=3, s=:dash, extra..., ytickfontcolor=c, xlims=[0,2pi])  # better?
        end
        data = getfield(nt, sy)

        data = vcat(data, data[1])  # wrap around
        # @show maximum(data)
        plot!(twinx(), theta, data; c, extra...,
            #yguide=string('\n'^cnt, sy), yguidefontcolor=c,
            ytickfontcolor=c, yticks, fill=0, fillalpha=0.1, xlims=[0,2pi], yguide=yguide2, yguidefontcolor=c)
        cnt += 1
    end
    γ = round(nt.gamma; sigdigits)
    s = round(nt.speed; sigdigits)
    info = round(nt.info; sigdigits)
    # title="γDᵣ=$γ ⇒ v/v₀=$s, i/Dᵣ=$info"
    title="i/Dᵣ=$info, v/v₀=$s"
    if hasproperty(nt, :tau) && nt.tau !== nothing
        tau = round(nt.tau; sigdigits)
        eta = round(nt.eta; sigdigits)
        title = title * "\nτDᵣ=$tau, η=$eta"
    end
    plot!(p; legend=false, title)
end

whole2plot(many::Vector; title="", size=(1200, 1000), lock::Bool=true, kw...) = plot(whole2plot.(many; lock)...; plot_title=title, size, kw...)


#####
##### Plotting -- complex strategies
#####

function whole2heatmap(nt::NamedTuple; both::Bool=true, color=cgrad([:lightgrey, :blue, :darkred]), fun::Function=identity, sigdigits=3)
    n = length(nt.prob)
    θs = angles2W(n)
    if !hasproperty(nt, :jump)
        if hasproperty(nt, :lambda)
            jump = zeros(n-1, n)
            jump .= (nt.lambda)' ./ (n-1)  # right? This preserves total rate... but some chance I want 2pi somewhere to make the numbers prettier.
            return whole2heatmap((; nt..., jump); both, color, fun, sigdigits)
        elseif hasproperty(nt, :zeta)
            jump = zeros(n-1, n)
            jump[n÷2, :] .= nt.zeta
            return whole2heatmap((; nt..., jump); both, color, fun, sigdigits)
        elseif hasproperty(nt, :mu)
            jump = steertojump(nt, 2).jump
            return whole2heatmap((; nt..., jump); both, color, fun, sigdigits)
        end
        error("you could make this translate, but not yet!")
    end
    # @show color
    p1 = heatmap(θs, θs[2:end], fun.(_full2Wjump(nt.jump)); yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yguide="change Δθ", yticks=([0, π, 2pi], ["0", "π", "2π"]),
        # color=:dense,
        color,
        # color=cgrad(:default),
        # color=cgrad([:lightgrey, :blue, :darkred]),
    )
    p2 = heatmap(θs, θs, fun.(jump2Wfold(nt.jump)); yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yguide="destination θ'", yticks=([0, π, 2pi], ["0", "π", "2π"]),
        color=:dense,
    )
    γ = round(nt.gamma; sigdigits)
    s = round(nt.speed; sigdigits)
    info = round(nt.info; sigdigits)
    if both
        plot(p1, p2; plot_title="γDᵣ=$γ ⇒ v/v₀=$s, i/Dᵣ=$info", size=(500,800), layout=(2,1))
    else
        # plot!(p1; title="γ=$γ\nv/v₀=$s, i/Dᵣ=$info\n\n")
        # plot!(p1; title="γDᵣ=$γ ⇒ v/v₀=$s, i/Dᵣ=$info")
        plot!(p1; title="i/Dᵣ=$info, v/v₀=$s")
    end
end

whole2heatmap(many::Vector; title="", size=(1200, 1000), color=cgrad([:lightgrey, :blue, :darkred]), sigdigits=3, kw...) = plot(whole2heatmap.(many; both=false, color, sigdigits)...;  plot_title=title, size, kw...)

function whole2diff(nt::NamedTuple; both::Bool=true, yguide="diffusion time")
    n = length(nt.prob)
    θs = angles2W(n)
    p1 = heatmap(θs, θs[2:end], _full2Wjump(nt.jump); yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yguide="change Δθ", yticks=([0, π, 2pi], ["0", "π", "2π"]),
        color=cgrad([:lightgrey, :blue, :darkred]),
    )
    # That's as above.
    # ts = size(nt.rates, 1)
    m = length(nt.taus)
    yticks = round.([nt.taus[1], nt.taus[1 + end÷2], nt.taus[end]]; sigdigits=2)
    yticks = (yticks, string.(yticks))
    p3 = heatmap(θs, nt.taus, nt.rates; yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yaxis=:log10,
        # yticks,
        yguide = "turning time (m=$m)",
        color=cgrad([:lightgrey, :blue, :darkred]),
    )

    γ = round(nt.gamma, sigdigits=2)
    Dt = nt.Dt isa Integer ? nt.Dt : round(nt.Dt, sigdigits=2)
    s = round(nt.speed, sigdigits=3)
    info = round(nt.info, sigdigits=3)
    eta = round(nt.eta, sigdigits=3)
    if both
        plot(p3, p1; plot_title="γ=$γ, Dₜ/Dᵣ=$Dt\nv/v₀=$s, i/Dᵣ=$info, η=$eta", size=(500,800), layout=(2,1))
    else
        plot!(p3; title="γ=$γ, Dₜ/Dᵣ=$Dt\nv/v₀=$s, i/Dᵣ=$info, η=$eta", yguide="time τDᵣ")
    end
end

whole2diff(many::Vector; title="", kw...) = plot(whole2diff.(many; both=false, kw...)...; size=(1200, 1000), plot_title=title)


function whole2box(nt::NamedTuple; both::Bool=true)
    n = length(nt.prob)
    θs = angles2W(n)
    p1 = heatmap(θs, θs[2:end], _full2Wjump(nt.jump); yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yguide="change Δθ", yticks=([0, π, 2pi], ["0", "π", "2π"]),
        color=cgrad([:lightgrey, :blue, :darkred]),
    )
    # That's as above.

    m = size(nt.rates, 1)
    p3 = heatmap(θs, θs[2:m+1], nt.rates; yflip=false,
        xguide="θ (n=$n)", xticks=([0, π, 2pi], ["0", "π", "2π"]),
        yguide = "max Δθ (m=$m)", yticks=([0, π], ["0", "π"]),
        color=cgrad([:lightgrey, :blue, :darkred]),
    )

    γ = round(nt.gamma, sigdigits=2)
    s = round(nt.speed, sigdigits=3)
    info = round(nt.info, sigdigits=3)
    if both
        plot(p3, p1; plot_title="γ=$γ\nv/v₀=$s, i/Dᵣ=$info", size=(500,800), layout=(2,1))
    else
        plot!(p3; title="γ=$γ, v/v₀=$s, i/Dᵣ=$info")
    end
end

whole2box(many::Vector; title="", kw...) = plot(whole2box.(many; both=false, kw...)...; size=(1200, 1000), plot_title=title)


#####
##### The end.
#####

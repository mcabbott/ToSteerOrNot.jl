export max2Whybrid, max2Watomic, contactfun, max2Wpair, max2Wfixpair
export max2Wcostumble, max2Wcosreverse, max2Wcosabs, max2Wsinsteer, max2Wcospair, max2Wcosflick, max3sinjump, max3deltajump, max3sindriftNEW


"""
    max2Whybrid(γ, n=20; steps=1000, mirror=true)

This is like `max2Wjump(; signed=false, mirror=true)`, except that
it finds Q by an inner Blahut-Arimoto-like iteration,
while using gradient-based `maximise` only for `lambda`.

With `log=false`, it runs and finds bad solutions.
With `log=true`, usually the solver gives up once the gradient gets a NaN,
which seems to happen as lambda becomes small, but not really tracked down!

```
gamma = 0.1
sol1 = max2Wjump(gamma, 40; signed=false, mirror=true)
sol2 = max2Whybrid(gamma, 20; log=true)
whole2heatmap([sol1, sol2]; size=(600,300))

sol2.qbar
sol2.contact
sol2.chi
sol2.lambda
sol2.calls  # outer iterations

sol2in = ToSteerOrNot.inner_BA_step(sol2.prob, sol2.lambda, sol2.gamma, 1000, sol2.mirror)

let sol2 = sol2
  n = length(sol2.lambda)
  x = angles2W(n)[2:end]
  p1 = plot(x, sol2.contact, lab="Psi(Δθ)", xguide="Δθ")
  plot!(p1, x, sol2.qbar, lab="qbar(Δθ)", fill=0, fillalpha=0.1)
  _wrap(x) = vcat(x, x[1])
  x2 = vcat(angles2W(n), 2pi)
  p2 = plot(x2, _wrap(sol2.chi), lab="χ(θ)", c=3, xguide="θ")
  plot!(p2, x2, _wrap(sol2.lambda ./ maximum(sol2.lambda)), lab="λ(θ), scaled", c=4, fill=0, fillalpha=0.1)
  plot!(p2, x2, _wrap((n/4) .* sol2.prob), lab="p(θ), scaled", c=:grey, fill=0, fillalpha=0.1)
  plot(p1, p2, size=(800,400), plot_title=string("gamma=", gamma, ", n=",n))
end
```

Variant allowing sign:
```
sol1s = max2Wjump(gamma, 40; signed=true, mirror=false)
sol2s = max2Whybrid(gamma, 40; mirror=false, log=true)
whole2heatmap([sol1s, sol2s]; size=(600,300))
```
"""
function max2Whybrid(gamma::Real, n::Int=20; init=ones, steps::Int=1000, mirror::Bool=true, kw...)
    lambda = init(n) ./ 10
    last_p = fill(1/n, n)
    call_count = Ref(0)
    maximise((; lambda); kw...) do nt
        call_count[] += 1

        # First we find qbar, using the p from last step:
        jump, qbar, chi, contact = inner_BA_step(last_p, nt.lambda, gamma, steps, mirror)
        # jump, qbar, chi, contact = inner_BA_step(last_p, ForwardDiff.value.(nt.lambda), gamma, steps, mirror)  # this is too much, doesn't solve for lambda at all

        @assert size(jump) == (n-1, n)
        # Then we compute FP and solve for p as usual:
        mat = master2Wjump(1, jump)

        if !(all(isfinite, mat))
            @warn "found Inf or NaN" call_count[]
            # display(ForwardDiff.value.(mat))
            println("jump = ")
            display(ForwardDiff.value.(jump))
            println("chi = ")
            display(ForwardDiff.value.(chi))
            println("qbar = ")
            display(ForwardDiff.value.(qbar))
        end

        prob = zero_eigvec_iter3(mat)
        last_p .= ForwardDiff.value.(prob)  # save for next iteration
        speed = cos2W(prob)
        info = info2jump(prob, jump)

        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, qbar, chi, contact, n, mirror)
    end
end

function inner_BA_step(prob::Vector, lambda::Vector, gamma::Real, steps::Int, mirror::Bool=true)
    n = length(lambda)

    ### Lagrange multiplier obeying 0 == χ''(θ) + cos(θ) + γ(λ(θ)-λ0)
    chi = inner_chi3(lambda, gamma)
    # chi = inner_chi3(ForwardDiff.value.(lambda), gamma)  # doesn't change much

    # if !(all(isfinite, chi))
    #     display(chi)
    #     display(lambda)
    #     error()
    # end

    ### rho, exp(χ-χ) matrix
    rho = prob .* lambda ./ dot(prob, lambda)
    _topmat = if mirror
        [(chi[mod1(μ+i, n)] + chi[mod1(μ-i, n)]) / (2*gamma) for μ in 1:n, i in 1:n-1]
    else
        [(chi[mod1(μ+i, n)]) / gamma for μ in 1:n, i in 1:n-1]
    end
    Emat = exp.(_topmat .- maximum(_topmat) .+ 10)
    @assert size(Emat) == (n, n-1)

    ### BA loop to find qbar
    qbar = fill!(similar(rho, n-1), 1/(n-1))  # initialise, should sum to 1
    # qbar[1] = 0
    # qbar[end] = 0
    # local zed, contact  # scope these outside the for loop
    # for _ in 1:steps
    #     zed = Emat * qbar
    #     contact = transpose(Emat) * (rho ./ zed)
    #     qbar = contact .* qbar
    # end
    # Faster version re-using same arrays:
    zed = similar(rho, n)
    contact = similar(rho, n-1)
    _rat = similar(zed)
    for ba_step in 1:steps
        mul!(zed, Emat, qbar)
        if any(iszero, zed)
            @info "got a zero in zed" ba_step
            println("zed = ")
            display(ForwardDiff.value.(zed))
            println("chi = ")
            display(ForwardDiff.value.(chi))
            println("qbar = ")
            display(ForwardDiff.value.(qbar))
            error()
        end
        _rat .= rho ./ (zed .+ 1e-12)
        mul!(contact, transpose(Emat), _rat)
        qbar .= contact .* qbar
    end

    if !(all(isfinite, qbar))
        @info "qbar contains Inf or NaN"
        println("chi = ")
        display(ForwardDiff.value.(chi))
        println("E = ")
        display(ForwardDiff.value.(Emat))
        @show extrema(Emat)
        println("qbar = ")
        display(ForwardDiff.value.(qbar))
        println("zed = ")
        display(ForwardDiff.value.(zed))
        error()
    end
    if any(isnan, ForwardDiff.partials.(qbar, 2))
        @info "qbar's gradient contains NaN"
        println("chi = ")
        display(ForwardDiff.value.(chi))
        println("qbar = ")
        display(ForwardDiff.value.(qbar))
    end

    ### jump is λ(Δθ|θ) = λ(θ) q_θ(Δθ)
    jump = transpose(lambda) .* transpose(Emat) .* qbar ./ transpose(zed)
    @assert size(jump) == (n-1, n)

    (; jump, qbar, chi, contact, Emat, zed)
end

function inner_chi1(lambda::Vector, gamma::Real)
    n = length(lambda)
    θs = angles2W(n)
    dθ = step(θs)
    @assert dθ ≈ 2pi/n

    # Lagrange multiplier obeying 0 == χ''(θ) + cos(θ) + γ(λ(θ)-λ0)
    # My first extremely crude implementation!
    lam0 = dθ/2pi * sum(lambda)
    rhs = @. -cos(θs) - gamma * (lambda - lam0)
    chiprime = cumsum(rhs).*dθ
    # chi = cumsum(chiprime .- mean(chiprime)).*dθ
    chi = cumsum(chiprime).*dθ
    @assert length(chi) == n
    # Then pick something to make it roughly periodic in θ?
    # c1 = chi[1]-chi[end]/θs[end]
    c1 = chi[1]-chi[end]/2pi
    chi .= chi .+ c1.*(θs)
end

# function inner_chi2(lambda::Vector, gamma::Real)
#     n = length(lambda)
#     θs = angles2W(n)
#     dθ = step(θs)
#     @assert dθ ≈ 2pi/n

#     lam0 = sum(lambda)/n  # dθ/2pi * sum(lambda), the same!
#     rhs = @. -cos(θs) - gamma * (lambda - lam0)

#     # This is what Jose does:
#     L2 = SymTridiagonal(fill(-2,n-1), fill(1,n-2))
#     chi = vcat(0, L2 \ (rhs[2:end] .* dθ^2))
# end

function inner_chi3(lambda::Vector, gamma::Real)
    n = length(lambda)
    θs = angles2W(n)
    dθ = step(θs)
    @assert dθ ≈ 2pi/n

    lam0 = sum(lambda)/n  # dθ/2pi * sum(lambda), the same!
    rhs = @. -cos(θs) - gamma * (lambda - lam0)

    # This seems more honest BC?
    L = collect(SymTridiagonal(fill(-2,n), fill(1,n-1)))
    L[1,n] = 1
    L[n,1] = 1

    tmp = (L \ rhs)

    chi = tmp .= (tmp .- mean(tmp)) .* dθ^2
end



"""
    max2Wpair(γ, n=20; mirror=true)

This is like `max2Wjump(; signed=false, mirror=true)`, except that
it optimises only `lambda` and `qbar`, using some expressions derived
by assuming we're at the optimum.

```
gamma = 0.005
gamma = 0.05
gamma = 0.5

sol1 = max2Wjump(gamma, 40; signed=false, mirror=true)  # brute force
sol2 = max2Wpair(gamma, 56; log=true)
whole2heatmap([sol1, sol2]; size=(600,300))

let
  dx1 = range(0,2pi,sol1.n+1)[2:end-1]
  plot(dx1, contactfun(sol1)[1], xguide="Δθ", lab="contact 1", c=1)
  plot!(dx1, contactfun(sol1)[2], lab="qbar(Δθ), 1", c=1, fill=0, fillalpha=0.1)

  dx2 = range(0,2pi,sol2.n+1)[2:end-1]
  plot!(dx2, contactfun(sol2)[1], xguide="Δθ", lab="contact 2", c=2)
  plot!(dx2, contactfun(sol2)[2], lab="qbar(Δθ), 2", c=2, fill=0, fillalpha=0.1)
end
```
"""
function max2Wpair(gamma::Real, n::Int=20; init=ones, signed::Bool=false, mirror::Bool=true, kw...)
    iseven(n) || error("haven't thought about odd n here")
    if signed
        lambda_raw = init(n) ./ 10
    else
        lambda_raw = init(n÷2 + 1) ./ 10
    end
    if mirror
        # negative is just a hack so that log=true only transforms the other parameter!
        minus_qbar = .-init(n÷2) ./ 10
    else
        minus_qbar = .-init(n-1) ./ 10
    end
    # last_p = fill(1/n, n)

    maximise((; lambda_raw, minus_qbar); kw...) do nt
        if signed
            lambda = nt.lambda_raw
        else
            lambda = vcat(nt.lambda_raw, @view nt.lambda_raw[reverse(2:end-1)])
        end
        if mirror
            half_qbar = @. clamp(-nt.minus_qbar, 0, Inf)
            qbar = vcat(half_qbar, @view half_qbar[reverse(1:end-1)])
        else
            qbar = @. clamp(-nt.minus_qbar, 0, Inf)
        end
        @assert length(lambda) == n
        @assert length(qbar) == n-1

        # Now reconstruct the full jump matrix... using the same pieces as for hybrid thing above.

        chi = inner_chi3(lambda, gamma)

        # rho = prob .* lambda ./ dot(prob, lambda)
        _topmat = if mirror
            [(chi[mod1(μ+i, n)] + chi[mod1(μ-i, n)]) / (2*gamma) for μ in 1:n, i in 1:n-1]
        else
            [(chi[mod1(μ+i, n)]) / gamma for μ in 1:n, i in 1:n-1]
        end
        Emat = exp.(_topmat .- maximum(_topmat) .+ 10)
        @assert size(Emat) == (n, n-1)
        zed = Emat * qbar

        # jump = transpose(lambda) .* transpose(Emat) .* qbar ./ transpose(zed)
        jump = transpose(lambda) .* transpose(Emat) .* qbar ./ (transpose(zed) .+ 1e-12)
        @assert size(jump) == (n-1, n)

        mat = master2Wjump(1, jump)

        if !(all(isfinite, mat))
            @info "mat contains Inf or NaN"
            println("chi = ")
            display(ForwardDiff.value.(chi))
            println("E = ")
            display(ForwardDiff.value.(Emat))
            @show extrema(Emat)
            println("qbar = ")
            display(ForwardDiff.value.(qbar))
            println("zed = ")
            display(ForwardDiff.value.(zed))
            error()
        end

        prob = zero_eigvec_iter3(mat)
        # last_p .= ForwardDiff.value.(prob)  # save for next iteration
        # last_p .= (1-alpha) .* last_p .+ alpha .* ForwardDiff.value.(prob)  # save for next iteration
        speed = cos2W(prob)
        info = info2jump(prob, jump)

        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, lambda, qbar, chi, n, mirror)
    end
end

"""
    max2Wcospair(γ, n=20; steps=1000, mirror=true, order=10)

Version of `max2Wpair(γ)` which uses a cosine basis for `λ(θ)`.
I tried quite a few things to try to make `qbar` behave well, don't know why.
"""
function max2Wcospair(gamma::Real, n::Int=20; order=10, kw...)
    iseven(n) || error("haven't thought about odd n here")
    # half_qbar = ones(n÷2) ./ 10
    half_qbar = ones(n-1) ./ 10

    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]

    maximise((; coscoeff, half_qbar); kw...) do nt
        lambda = clamp.(cosmat * nt.coscoeff, 0, Inf)

        # _qbar = vcat(nt.half_qbar, @view nt.half_qbar[reverse(1:end-1)])
        _qbar = @views nt.half_qbar .+ nt.half_qbar[reverse(1:end)]
        qbar = _qbar ./ sum(_qbar)  # IDK why I didn't do this above, although it ought to cancel

        @assert length(lambda) == n
        @assert length(qbar) == n-1

        chi = inner_chi3(lambda, gamma)

        _topmat = [(chi[mod1(μ+i, n)] + chi[mod1(μ-i, n)]) / (2*gamma) for μ in 1:n, i in 1:n-1]

        Emat = exp.(_topmat .- maximum(_topmat) .+ 10)
        @assert size(Emat) == (n, n-1)
        zed = Emat * qbar

        jump = transpose(lambda) .* transpose(Emat) .* qbar ./ (transpose(zed) .+ 1e-12)
        @assert size(jump) == (n-1, n)

        mat = master2Wjump(1, jump)

        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2jump(prob, jump)

        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, lambda, qbar, chi, n, mirror=true, signed=false)
    end
end

function max2Wfixpair(gamma::Real, n::Int=20; order=10, kw...)
    iseven(n) || error("haven't thought about odd n here")

    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]

    atoms=4
    n < 16 && error("can't do n<16 with 4 atoms")
    qbarcoeff = ones(atoms) ./ atoms

    maximise((; coscoeff, qbarcoeff); kw...) do nt
        lambda = clamp.(cosmat * nt.coscoeff, 0, Inf)

        qbar = zeros(eltype(lambda), n-1)
        qbar[n÷2] = nt.qbarcoeff[1]
        qbar[n÷4] = qbar[n - n÷4] = nt.qbarcoeff[2]
        qbar[n÷8] = qbar[n - n÷8] = nt.qbarcoeff[3]
        qbar[n÷16] = qbar[n - n÷16] = nt.qbarcoeff[4]

        @assert length(lambda) == n
        @assert length(qbar) == n-1

        chi = inner_chi3(lambda, gamma)

        _topmat = [(chi[mod1(μ+i, n)] + chi[mod1(μ-i, n)]) / (2*gamma) for μ in 1:n, i in 1:n-1]

        Emat = exp.(_topmat .- maximum(_topmat) .+ 10)
        @assert size(Emat) == (n, n-1)
        zed = Emat * qbar

        jump = transpose(lambda) .* transpose(Emat) .* qbar ./ (transpose(zed) .+ 1e-12)
        @assert size(jump) == (n-1, n)

        mat = master2Wjump(1, jump)

        prob = zero_eigvec_iter3(mat)
        speed = cos2W(prob)
        info = info2jump(prob, jump)

        (; objective = speed - gamma * info, speed, info, gamma, prob, jump, lambda, qbar, chi, n, mirror=true, signed=false)
    end
end

"""
    contactfun(sol)

Returns the thing formally known as `Γ(Δθ)`, given the output of `max2Wjump`.
Better version of `newGamma`.

```
gamma = 0.5

gamma = 0.005
gamma = 0.05
sol = max2Wjump(gamma, 40; signed=false, mirror=true)
con = contactfun(sol)

let p1 = whole2heatmap(sol, both=false)
  dx = range(0,2pi,sol.n+1)[2:end-1]
  p2 = plot(dx, con[1], xguide="Δθ", lab="contact fn")
  plot!(p2, dx, con[2], lab="qbar(Δθ)")
  x = range(0,2pi,sol.n+1)[1:end-1]
  p3 = plot(x, con[3], xguide="θ", lab="λ(θ)")
  plot!(p3, x, con[4], lab="χ(θ)")
  plot(p1, p2, size=(800,400))
end

```
"""
function contactfun(nt::NamedTuple)
    n = length(nt.prob)
    @assert size(nt.jump) == (n-1, n)
    θs = angles2W(n)
    dθ = step(θs)
    @assert dθ ≈ 2pi/n

    p = nt.prob  # p(θ) dθ
    @assert sum(p) ≈ 1
    lambda = vec(sum(nt.jump; dims=1))  # λ(θ)
    # jump is λ(Δθ|θ) = λ(θ) q(Δθ|θ), work out this q matrix, Jose calls it q_θ(Δθ)
    qmat = nt.jump ./ (transpose(lambda) .+ 1e-10)  # 1e-10

    rho = p .* lambda ./ dot(p, lambda)  # ρ(θ) dθ
    @assert length(rho) == n

    qlam = qmat * rho  # q^λ(Δθ) = int dθ ρ(θ) q(Δθ|θ)  sometimes qbar in notes
    @assert length(qlam) == n-1

    chi = inner_chi3(lambda, nt.gamma)
    _topmat = if nt.mirror
        [(chi[mod1(μ+i, n)] + chi[mod1(μ-i, n)]) / (2*nt.gamma) for μ in 1:n, i in 1:n-1]
    else
        [(chi[mod1(μ+i, n)]) / nt.gamma for μ in 1:n, i in 1:n-1]
    end
    Emat = exp.(_topmat .- maximum(_topmat) .+ 10)
    zed = Emat * qlam
    Psi = transpose(Emat) * (rho ./ zed)

    (; Psi, qlam, lambda, chi)
end


export newGamma



#=

Using a Fourier basis for numerical solutions,
instead of brute-force solving... works very well!
Should use for final plots.

=#


"""
    max2Wcostumble(γ, n=72; tau=nothing)

This is like `max2Wtumble`, but works in a basis of cosine functions, instead of directly `λ(θ)`.
Better & faster!

```
[max2Wcostumble(g; order=10) for g in onethree(0.01, 3)] |> reverse |> whole2plot;
plot!(size=(800,600))
```
"""
function max2Wcostumble(gamma::Real, n::Int=72; tau=nothing, order=10, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]
    maximise((; coscoeff); kw...) do nt
        lambda = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master2Wtrio(1, lambda, nothing, nothing)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        eta = eta2(prob, lambda, tau)
        speed = cos2W(prob) * eta
        info = info2rate(prob, lambda)
        (; objective = speed - gamma * info, speed, info, gamma, prob, lambda, tau, eta)
    end
end


"""
    max2Wcosreverse(γ, n=72; tau=nothing)

This is like `max2Wreverse`, but works in a basis of cosine functions.
Better & faster!

```
[max2Wcosreverse(g; order=10) for g in onethree(0.01, 3)] |> reverse |> whole2plot;
plot!(size=(800,600))
```
"""
function max2Wcosreverse(gamma::Real, n::Int=72; tau=nothing, order=10, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]
    maximise((; coscoeff); kw...) do nt
        zeta = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master2Wreverse(1, zeta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        eta = eta2(prob, zeta, tau)
        speed = cos2W(prob) * eta
        info = info2rate(prob, zeta)
        (; objective = speed - gamma * info, speed, info, gamma, prob, zeta, tau, eta)
    end
end


"""
    max2Wcosflick(γ, n=72; tau=nothing)

This is like `max2Wflick`, but works in a basis of cosine functions.
Better & faster!

```
[max2Wcosreverse(g; order=10) for g in onethree(0.01, 3)] |> reverse |> whole2plot;
plot!(size=(800,600))
```
"""
function max2Wcosflick(gamma::Real, n::Int=72; tau=nothing, order=10, delta=div(n, 4), kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]
    maximise((; coscoeff); kw...) do nt
        kappa = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master2Wtrio(1, nothing, kappa, nothing; delta)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        eta = eta2(prob, kappa, tau)
        speed = cos2W(prob) * eta
        info = info2rate(prob, kappa)
        (; objective = speed - gamma * info, speed, info, gamma, prob, kappa, tau, eta)
    end
end

"""
    max2Wcosabs(γ, n=72; tau=nothing)

This is like `max2Wdrift(γ, n; signed=false)`, but works in a basis of cosine functions.
Better & faster!

```
[max2Wcosabs(g; order=30) for g in onethree(0.01, 3)] |> reverse |> whole2plot;
plot!(size=(800,600))
```
"""
function max2Wcosabs(gamma::Real, n::Int=72; order=30, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [cos(k*θ) for θ in θs, k in 0:order-1]
    maximise((; coscoeff, Dc=1.0); kw...) do nt
        mu = cosmat * nt.coscoeff
        mat = master2Wtrio(1+nt.Dc, nothing, nothing, mu)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        speed = cos2W(prob)
        speedx = sin2W(prob)
        info = info2drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, gamma, prob, mu, signed=false)
    end
end

"""
    max2Wsinsteer(γ, n=72; tau=nothing)

This is like `max2Wdrift(γ, n; signed=true)`, but works in a basis of sin functions.
Better & faster!
```
[max2Wsinsteer(g; order=5) for g in onethree(0.01, 3)] |> reverse |> whole2plot;
plot!(size=(800,600))
```
"""
function max2Wsinsteer(gamma::Real, n::Int=72; order=5, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles2W(n)
    cosmat = [sin(k*θ) for θ in θs, k in 1:order]
    maximise((; coscoeff, Dc=1.0); kw...) do nt
        mu = cosmat * nt.coscoeff
        mat = master2Wtrio(1+nt.Dc, nothing, nothing, mu)
        # prob = zero_eigvec(mat)
        prob = zero_eigvec_iter3(mat; shift=1e-10)
        speed = cos2W(prob)
        speedx = sin2W(prob)
        info = info2drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, gamma, prob, mu, signed=true)
    end
end






export max3sintumble, max3sinflick, max3sinreverse, max3sindrift

"""
    max3sintumble(γ, n=36)

Version of `max3tumble` using sine basis.
"""
function max3sintumble(gamma::Real, n::Int=36; order=5, tau=nothing, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]
    res = maximise((; coscoeff); kw...) do nt
        lambda = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master3d(1, lambda, nothing)
        prob = zero_eigvec_iter3(mat, norm3d)
        eta = eta3(prob, lambda, tau)
        speed = cos3d(prob) * eta
        info = info3d(prob, lambda)
        (; objective = speed - gamma * info, speed, info, prob, lambda, gamma, eta)
    end
end


"""
    max3sinflick(γ, n=36)

Version of `max3flick` using sine basis.
"""
function max3sinflick(gamma::Real, n::Int=36; order=5, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]
    res = maximise((; coscoeff); kw...) do nt
        kappa = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master3d(1, nothing, kappa)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3d(prob, kappa)
        (; objective = speed - gamma * info, speed, info, prob, kappa, gamma)
    end
end

"""
    max3sinreverse(γ, n=36)

Version of `max3reverse` using sine basis.
"""
function max3sinreverse(gamma::Real, n::Int=36; order=5, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = 0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]
    res = maximise((; coscoeff); kw...) do nt
        zeta = clamp.(cosmat * nt.coscoeff, 0, Inf)
        mat = master3d(1, nothing, nothing, zeta)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3d(prob, zeta)
        (; objective = speed - gamma * info, speed, info, prob, zeta, gamma)
    end
end

"""
    max3sindrift(γ, n=36)

Version of `max3drift` using sine basis.
"""
function max3sindrift(gamma::Real, n::Int=36; order=5, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = -0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]
    res = maximise((; coscoeff, Dc=1.0); kw...) do nt
        mu = cosmat * nt.coscoeff
        mat = master3drift(1 + nt.Dc, mu)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, prob, mu, gamma, signed=true)
    end
end

function max3sindriftNEW(gamma::Real, n::Int=36; order=5, kw...)
    coscoeff = zeros(order)
    coscoeff[1] = -0.1
    θs = angles3(n)
    sinθs = sin.(θs)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]
    res = maximise((; coscoeff, Dc=1.0); kw...) do nt
        mu = cosmat * nt.coscoeff
        mat = master3drift(1, nt.Dc, mu, sinθs)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3drift(prob, nt.Dc, mu)
        (; objective = speed - gamma * info, speed, info, prob, mu, gamma, signed=true)
    end
end

"""
    max3sinjump(γ, n=20, m=n)

Version of `max3jump` using sine basis.

```
sol_1s = max3sinjump(1.0, 16; iterations=1000)
jump3plot(sol_1s)

sol_2s = max3sinjump(0.1, 16; iterations=1000)
jump3plot(sol_2s)  # still reverse

sol_3s = max3sinjump(0.01, 16; iterations=1000)
jump3plot(sol_3s)  # reverse + flick

# sol_4a = max3sinjump(0.001, 16; iterations=3000)  # 6 seconds
# sol_4b = max3sinjump(0.001, 20; iterations=3000)  # 5 minutes!
sol_4s = max3sinjump(0.001, 16, 32; iterations=1000, order=10)
jump3plot(sol_4s, fun=sqrt)  # maybe more?

plot(jump3plot(sol_1s), jump3plot(sol_2s), jump3plot(sol_3s), jump3plot(sol_4s); size=(800, 700))
```
"""
function max3sinjump(gamma::Real, n::Int=20, m::Int=n; order=5, kw...)
    # jump = init(m, n)/10
    coscoeff = ones(m, order) ./ 10
    coscoeff[1] = -0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]'

    deltas = range(0, pi, length=m+1)[2:end]

    sinθs, cosθs = eachcol(stack(sincos, θs; dims=1))
    sinΔs, cosΔs = eachcol(stack(sincos, deltas; dims=1))
    ψs = range(0, 2pi, length=4n+1)[1:end-1]
    sinψs, cosψs = eachcol(stack(sincos, ψs; dims=1))

    res = maximise((; coscoeff); kw...) do nt
        # jump = clamp.(nt.coscoeff * cosmat, 0, Inf)
        jump = exp.(nt.coscoeff * cosmat)
        mat = master3jump(1, jump, deltas, sinθs, cosθs, sinψs, cosψs, sinΔs, cosΔs)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3jump(prob, jump)
        (; objective = speed - gamma * info, speed, info, prob, gamma, deltas, jump)
    end
end


"""
    max3deltajump(γ, n=20, m=4)

Version of `max3jump` using sine basis, and also adjusting deltas.
Doesn't work very well, many cases run into `Inf` and I'm not sure why yet.

```
sol_2d = max3deltajump(0.1, 32; iterations=2000)
jump3plot(sol_2d)  # reverse

sol_3d = max3deltajump(0.01, 30; iterations=1000)
jump3plot(sol_3d)  # reverse + flick

sol_4d = max3deltajump(0.001, 16, 6; iterations=1000, order=10)
jump3plot(sol_4d, fun=sqrt)
```
"""
function max3deltajump(gamma::Real, n::Int=20, m::Int=4; order=5, kw...)

    coscoeff = ones(m, order) ./ 10
    coscoeff[1] = -0.1
    θs = angles3(n)
    cosmat = [k==0 ? 1.0 : sin(k*θ/2) for θ in θs, k in 0:order-1]'

    raw_deltas = collect(range(0, pi, length=m+1)[2:end])

    sinθs, cosθs = eachcol(stack(sincos, θs; dims=1))
    ψs = range(0, 2pi, length=4n+1)[1:end-1]
    sinψs, cosψs = eachcol(stack(sincos, ψs; dims=1))

    res = maximise((; coscoeff, raw_deltas); kw...) do nt
        deltas = clamp.(nt.raw_deltas, 1e-5, pi-1e-5)
        sinΔs, cosΔs = eachcol(stack(sincos, deltas; dims=1))

        jump = exp.(nt.coscoeff * cosmat)
        mat = master3jump(1, jump, deltas, sinθs, cosθs, sinψs, cosψs, sinΔs, cosΔs)
        if any(!isfinite, mat)
            @warn "about to fail!" deltas
            display(jump)
        end
        # prob = zero_eigvec(mat, norm3d)
        prob = zero_eigvec_iter3(mat, norm3d)
        speed = cos3d(prob)
        info = info3jump(prob, jump)
        (; objective = speed - gamma * info, speed, info, prob, gamma, deltas, jump)
    end
end

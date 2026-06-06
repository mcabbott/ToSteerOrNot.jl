export info2rate, info2drift, info2jump, plogp, plogpoverq, COLOURS2D, colour2


info2rate(prob::AbstractVector, rate::Nothing) = false
function info2rate(prob::AbstractVector, rate::AbstractVector)
    expect = dot(prob, rate)
    sum(prob .* plogpoverq.(rate, expect))
end

function info2drift(prob::Vector, Diff::Real, mu::Vector, reflect::Bool=false)
    if reflect  # the other half has opposite sign, hence mean is zero
        varmu = sum(prob .* mu.^2)
    else  # we have the whole circle, or something equivalent
        # mumu = mean(mu) # unweighted mean, WRONG forever?
        mumu = dot(prob, mu)  # weighted mean of mu using p
        varmu = sum(prob .* (mu .- mumu).^2)
    end
    varmu/(4*Diff)
end

function info2jump(prob::Vector, jumprate::Matrix)
    pdelta = jumprate * prob
    # jumprate * prob' * log(jumprate / pdelta)
    sum(@. prob' * plogpoverq(jumprate, pdelta))
end

plogp(p::Real) = iszero(p) ? zero(p) : p * log(p)
plogpoverq(p::Real, q::Real) = (p*q<=0) ? zero(p*q) : p * log(p/q)


COLOURS2D = Dict(:lambda => :blue,  # tumble
    :mu => :red,  # drift / steering
    :sign => :darkred,  # Matt's sign-only thing
    :abs => :darkorange,  # for unsigned steering?
    :kappa => :purple,  # :darkcyan, # flick
    :zeta => :green,  # reverse
)
function colour2(nt::NamedTuple)
    for sy in keys(COLOURS2D)
        hasproperty(nt, sy) || continue
        if sy === :mu && nt.signed === false
            return COLOURS2D[:abs]
        end
        return COLOURS2D[sy]
    end
end

function _find_ylim(nt)
    ymax = 0.0
    ymin = 0.0
    for sy in keys(COLOURS2D)
        hasproperty(nt, sy) || continue
        ymax = max(ymax, maximum(getfield(nt, sy)))
        ymin = min(ymin, minimum(getfield(nt, sy)))
    end
    (; ylims=[1.1*ymin, 1.1*ymax])
end

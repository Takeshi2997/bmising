module Func
    include("./setup.jl")
    using .Const, LinearAlgebra

    function updateσ(zh)

        σ = -ones(Float32, Const.n)
        prob = 1 ./ (1 .+ exp.(-2.0 * zh))
        pup = rand(Float32, Const.n)
        for ix in 1:Const.n
            if pup[ix] < prob[ix]
                σ[ix] = 1.0
            end
        end
        return σ
    end

    function updateh(zσ)

        h = -ones(Float32, Const.n)
        prob = 1 ./ (1 .+ exp.(-2.0 * zσ))
        pup = rand(Float32, Const.n)
        for ix in 1:Const.n
            if pup[ix] < prob[ix]
                h[ix] = 1.0
            end
        end
        return h
    end
end

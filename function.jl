module Func
    include("./setup.jl")
    using .Const, LinearAlgebra

    function updateσ(h, w, α)

        σ = -ones(Float32, Const.n)
        z = (w * h .- 10.0^(-7)) .* α
        prob = exp.(z) ./ 2.0 ./ cosh.(z)
        pup = rand(Float32, Const.n)
        for ix in 1:Const.n
            if pup[ix] < prob[ix]
                σ[ix] = 1.0
            end
        end
        return σ
    end

    function updateh(σ, w, α)

        h = -ones(Float32, Const.n)
        z = transpose(w) * σ .* α
        prob = exp.(z) ./ 2.0 ./ cosh.(z)
        pup = rand(Float32, Const.n)
        for ix in 1:Const.n
            if pup[ix] < prob[ix]
                h[ix] = 1.0
            end
        end
        return h
    end
end

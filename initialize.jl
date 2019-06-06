module Init
    include("./setup.jl")
    using .Const, LinearAlgebra

    function σ()

        σ = ones(Float64, Const.n)
        return σ
    end

    function α(β)

        return acosh(exp(2.0 * β)) ./ 2.0
    end

    function w()

        weight = Array(Diagonal(fill(1.0, (Const.n, Const.n))))
        for ix in 1:Const.n
            index = findfirst(isequal(ix), Const.x)
            i = index[1]
            j = index[2]
            iu = mod(i + 1 - 1, Const.l) + 1
            ju = mod(j + 1 - 1, Const.l) + 1
            weight[Const.x[iu, j], ix] = 1.0
            weight[Const.x[i, ju], ix] = 1.0
        end
        return weight
    end
end

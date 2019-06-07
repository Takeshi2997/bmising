module Init
    include("./setup.jl")
    using .Const, LinearAlgebra

    function σ()

        σ = ones(Float64, Const.n)
        return σ
    end

    function α(β)

        return acosh(exp(2.0 * β)) / 2.0
    end

    function w()

        weight = Array(Diagonal(fill(1.0, (Const.n, Const.n))))
        for i in 1:Const.l - 1
            for j in 1:Const.l
                weight[Const.x[i + 1, j], Const.x[i, j]] = 1.0
            end
        end
        for j in 1:Const.l - 1
            for i in 1:Const.l
                weight[Const.x[i, j + 1], Const.x[i, j]] = 1.0
            end
        end
        return weight
    end
end

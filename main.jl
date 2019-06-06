include("./setup.jl")
include("./initialize.jl")
include("./function.jl")
using .Const, .Init, .Func, LinearAlgebra

w = Init.w()

f = open("magnetization.txt", "w")
g = open("logmagnetization.txt", "w")
for it in 1:500
    t = 0.01 * it
    β = 1 / t
    m = 0.0
    n = 0.0
    σ = Init.σ()
    α = Init.α(β)
    for i in 1:Const.iters_num
        h = Func.updateh(σ, w, α)
        σ = Func.updateσ(h, w, α)
        if i > Const.burnintime
            if i % Const.sample_interval == 0
                n += 1.0
                m += sum(σ) / Const.n
            end
        end
    end

    m /= n

#    if 2.269185314 - t > 0
#        write(g, string(log(2.269185314 - t)))
#        write(g, "\t")
#        write(g, string(log(m)))
#        write(g, "\n")
#    end
    write(f, string(t))
    write(f, "\t")
    write(f, string(m))
    write(f, "\n")
end
close(f)
close(g)

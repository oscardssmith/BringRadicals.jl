using BringRadical
using Test

@testset "BringRadical.jl" begin
    f(x) = x^5 + x
    for T in (Float32, Float64, BigFloat)
        for _ in 2^8
            for x in (rand(T), 20*rand(T), 100*rand(T))
                @test abs(f(bringrad(x)) - x) < eps(x)
                @test abs(f(bringrad(-x)) + x) < eps(x)
            end
        end
    end
end

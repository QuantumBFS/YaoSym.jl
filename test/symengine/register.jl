using Test, YaoSym, YaoBase

@testset "constructors" begin
    @test ket"011" isa SymReg
    @test (ket"011")' isa AdjointSymReg
end

@testset "partial_tr" begin
    @vars α β
    @test α * ket"101" + β * ket"111" |> partial_tr(2:3) == (α + β) * ket"1"
end

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    B = rand(10, 3)
    @compile_workload begin
        mva_eigen(B; dim = 1)
        mva(B; dim = 1)
        mva(B; dim = 1, field = :E)
        convection_efield(rand(10, 3), rand(10, 3); dim = 1)
    end
end

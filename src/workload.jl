using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    B = rand(10, 3)
    @compile_workload begin
        mva(B)
        mva(B; field = :E)
        convection_efield(rand(10, 3), rand(10, 3))
    end
end

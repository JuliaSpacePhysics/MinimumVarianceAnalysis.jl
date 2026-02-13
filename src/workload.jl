using PrecompileTools: @setup_workload, @compile_workload
using Unitful: nT

@setup_workload begin
    B = rand(10, 3) .* nT
    @compile_workload begin
        mva(B)
        mvae(B)
        convection_efield(rand(10, 3), rand(10, 3))
    end
end
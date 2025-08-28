using PrecompileTools: @setup_workload, @compile_workload
using Unitful: nT

@setup_workload begin
    B = rand(10, 3) .* nT
    @compile_workload begin
        mva(B)
    end
end
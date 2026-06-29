include(joinpath(@__DIR__, "binconv_utils.jl"))

const SRC = joinpath(@__DIR__, "..", "..", "Marked_edges")

# fig4_mew_cuts.csv / samples1.csv / cuts.csv hold integer cut-edge counts
# (small range, ~20-150) -> Int16. mm.csv holds Float64 percentages -> Float32
# (parse directly to Float32, no intermediate Float64 array).
convert_csv_to_binary(joinpath(SRC, "fig4_mew_cuts.csv"), joinpath(@__DIR__, "fig4_mew_cuts.bin"), Int16, 10)
convert_csv_to_binary(joinpath(SRC, "samples1.csv"), joinpath(@__DIR__, "fig4_baseline_samples.bin"), Int16, 1000)
convert_csv_to_binary(joinpath(SRC, "cuts.csv"), joinpath(@__DIR__, "fig5_6_cuts.bin"), Int16, 10)
convert_csv_to_binary(joinpath(SRC, "mm.csv"), joinpath(@__DIR__, "fig5_6_mm.bin"), Float32, 10)

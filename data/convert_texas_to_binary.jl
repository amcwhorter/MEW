include(joinpath(@__DIR__, "binconv_utils.jl"))

const SRC = joinpath(@__DIR__, "..", "..", "Texas_Optimization")

# Figure 8: 10 chains, intended at 4,000,000 steps each, but one chain
# (row 5) stopped at ~2,572,000 -- ragged rows are handled by binconv_utils.
# cuts.csv = cut-edge counts (Int16 is plenty, target range is ~3340-3420)
# -> Int16. mm.csv = mean-median score (range roughly -0.01 to 0.01) -> Float32.
convert_csv_to_binary(joinpath(SRC, "cuts.csv"), joinpath(@__DIR__, "fig8_tx_cuts.bin"), Int16, 10)
convert_csv_to_binary(joinpath(SRC, "mm.csv"), joinpath(@__DIR__, "fig8_tx_mm.bin"), Float32, 10)

# Regenerates Figures 3-8 from the compact data in data/. Each section is
# self-contained and frees its data before moving on (this machine has 8GB
# RAM). Run with: julia make_all_figures.jl
#
# Confirmed by eyeballing against the paper / plot_37.svg / direct user
# correction (see git history of this file for the back-and-forth):
#   - Figure 3 source: a re-run MEW chain (10M steps, 98.6% coverage), since
#     the original bludef.jl had a latent bug -- see fig3_enumerated/
#     fig3_sampled below and bludef_cheshire_save.jl in the laptop-side
#     Marked_edges/ folder for how it was regenerated.
#   - Figure 7(a) source: fig7_meangrid_ks_results4.csv.
#   - Figure 7(b) source: ks_var_10.csv (copied in as fig7_weightgrid_ks.csv).
#   - Figure 5(b) = ks_pairwise (mm) in blue, ks_to_target (originally "_a")
#     in red. Figure 6(b) = ks_pairwise (cuts) in blue, ks_to_target
#     (originally "_b") in red.

using Plots, Statistics, DelimitedFiles, JSON, Graphs, Shapefile, DataFrames
include(joinpath(@__DIR__, "data", "load_data.jl"))
include(joinpath(@__DIR__, "beano2.3.jl"))

const DATA = joinpath(@__DIR__, "data")
const OUT = joinpath(@__DIR__, "figures")
mkpath(OUT)

gaussian(x, beta, mu) = sqrt(beta / pi) * exp(-beta * (x - mu)^2)

println("=== Figure 3: Cheshire County (enumerated vs. sampled vs. unsampled) ===")
let
    data = JSON.parsefile(joinpath(@__DIR__, "NH", "NH_dual_graph_stripped.json"))
    nodes = data["nodes"]
    links = data["links"]
    cheshire_indices = collect(findall(node -> node["COUNTY"] == "Cheshire", nodes))
    g = Graphs.SimpleGraph(length(cheshire_indices))
    for link in links
        u, v = link["source"], link["target"]
        a = findfirst(x -> x == u + 1, cheshire_indices)
        b = findfirst(x -> x == v + 1, cheshire_indices)
        if !(isnothing(a) || isnothing(b))
            add_edge!(g, simple_edge(a, b))
        end
    end

    enum_rows = load_binary_gz(joinpath(DATA, "fig3_enumerated.bin.gz"), Int8)
    sampled_rows = load_binary_gz(joinpath(DATA, "fig3_sampled.bin.gz"), Int8)
    sampled_set = Set(sampled_rows)
    unsampled_rows = filter(r -> !(r in sampled_set), enum_rows)

    cuts_of(rows) = [length(cut_edges(r, g)) for r in rows]
    enum_cuts = cuts_of(enum_rows)
    sampled_cuts = cuts_of(sampled_rows)
    unsampled_cuts = cuts_of(unsampled_rows)

    pct = round(100 * length(sampled_rows) / length(enum_rows), digits=1)
    println("  coverage: $(pct)%")
    bin_edges = 9.5:1:37.5
    p = plot(xlabel="Cut Edges", ylabel="Count")
    histogram!(p, enum_cuts, bins=bin_edges, label="Enumerated")
    histogram!(p, sampled_cuts, bins=bin_edges, label="Sampled")
    histogram!(p, unsampled_cuts, bins=bin_edges, label="Unsampled")

    savefig(p, joinpath(OUT, "figure3.png"))
    println("  wrote figures/figure3.png")
end
GC.gc()

# Splits a (nrows x 2k) matrix from the old incremental-KS scripts into
# (iterations, ks_values), each (nrows x k), and returns the mean KS curve.
function mean_ks_curve(path::String)
    data = readdlm(path, ',')
    nrows, ncols = size(data)
    k = ncols ÷ 2
    iters = data[1, 1:k]
    vals = data[:, k+1:end]
    return iters, vec(mean(vals, dims=1))
end

println("=== Figure 4: spanning tree distribution vs. baseline ===")
let
    mew_rows = load_binary_gz(joinpath(DATA, "fig4_mew_cuts.bin.gz"), Int16)
    baseline_rows = load_binary_gz(joinpath(DATA, "fig4_baseline_samples.bin.gz"), Int16)

    p1 = plot(xlabel="Cut Edges", ylabel="Density")
    for (i, row) in enumerate(mew_rows)
        histogram!(p1, row, bins=19.5:1:140.5, normalize=true, fill=true, fillcolor=:lightblue,
                   linecolor=:blue, alpha=0.3, linealpha=0.5, label=(i == 1 ? "MEW (10 Runs)" : ""))
    end
    baseline_flat = reduce(vcat, baseline_rows)
    histogram!(p1, baseline_flat, bins=19.5:1:140.5, normalize=true, fill=nothing,
               linecolor=:red, linealpha=0.5, label="Tree Distribution")

    mew_rows = nothing; baseline_rows = nothing; baseline_flat = nothing; GC.gc()

    pair_iters, pair_ks = mean_ks_curve(joinpath(DATA, "fig4_ks_pairwise.csv"))
    base_iters, base_ks = mean_ks_curve(joinpath(DATA, "fig4_ks_to_baseline.csv"))
    p2 = plot(xlabel="Iteration", ylabel="KS Distance", xaxis=:log10, yaxis=:log10)
    plot!(p2, pair_iters, pair_ks, label="Average Pairwise KS Distance")
    plot!(p2, base_iters, base_ks, label="Average KS Distance to Tree Distribution")

    savefig(plot(p1, p2, layout=(1, 2), size=(1000, 400)), joinpath(OUT, "figure4.png"))
    println("  wrote figures/figure4.png")
end
GC.gc()

println("=== Figure 5: competitiveness ===")
let
    mm_rows = load_binary_gz(joinpath(DATA, "fig5_mm.bin.gz"), Float32)
    p1 = plot(xlabel="Democratic Percentage", ylabel="Density")
    for (i, row) in enumerate(mm_rows)
        histogram!(p1, row, bins=200, normalize=true, fill=true, fillcolor=:lightblue,
                   linecolor=:blue, alpha=0.3, linealpha=0.5, label=(i == 1 ? "MEW (10 Runs)" : ""))
    end
    plot!(p1, 49:0.01:51, gaussian.(49:0.01:51, 10, 50), color=:red, label="Target")
    mm_rows = nothing; GC.gc()

    pair_iters, pair_ks = mean_ks_curve(joinpath(DATA, "fig5_ks_pairwise.csv"))
    target_iters, target_ks = mean_ks_curve(joinpath(DATA, "fig5_ks_to_target.csv"))
    p2 = plot(xlabel="Iteration", ylabel="KS Distance", xaxis=:log10, yaxis=:log10)
    plot!(p2, pair_iters, pair_ks, label="Average Pairwise KS Distance")
    plot!(p2, target_iters, target_ks, label="Average KS Distance to Target")

    savefig(plot(p1, p2, layout=(1, 2), size=(1000, 400)), joinpath(OUT, "figure5.png"))
    println("  wrote figures/figure5.png")
end
GC.gc()

println("=== Figure 6: compactness ===")
let
    cuts_rows = load_binary_gz(joinpath(DATA, "fig6_cuts.bin.gz"), Int16)
    p1 = plot(xlabel="Cut Edges", ylabel="Density")
    for (i, row) in enumerate(cuts_rows)
        histogram!(p1, row, bins=64.5:1:85.5, normalize=true, fill=true, fillcolor=:lightblue,
                   linecolor=:blue, alpha=0.3, linealpha=0.5, label=(i == 1 ? "MEW (10 Runs)" : ""))
    end
    plot!(p1, 65:0.1:85, gaussian.(65:0.1:85, 0.1, 72), color=:red, label="Target")
    cuts_rows = nothing; GC.gc()

    pair_iters, pair_ks = mean_ks_curve(joinpath(DATA, "fig6_ks_pairwise.csv"))
    target_iters, target_ks = mean_ks_curve(joinpath(DATA, "fig6_ks_to_target.csv"))
    p2 = plot(xlabel="Iteration", ylabel="KS Distance", xaxis=:log10, yaxis=:log10)
    plot!(p2, pair_iters, pair_ks, label="Average Pairwise KS Distance")
    plot!(p2, target_iters, target_ks, label="Average KS Distance to Target")

    savefig(plot(p1, p2, layout=(1, 2), size=(1000, 400)), joinpath(OUT, "figure6.png"))
    println("  wrote figures/figure6.png")
end
GC.gc()

println("=== Figure 7: NH multivariate heatmaps ===")
let
    results4 = readdlm(joinpath(DATA, "fig7_meangrid_ks_results4.csv"), ',')
    demvals = 40:2:60   # best-effort, see caveats at top of file
    cutvals = 10:10:200
    p1 = heatmap(cutvals, demvals, results4, xlabel="Cut Edges", ylabel="Democratic Vote %")

    # Confirmed source: ks_var_10.csv. Axis order/clims match push_reader.jl's
    # actual heatmap(bet_values, alph_values, data, ..., clims=(0.1,1)) call.
    weightgrid = readdlm(joinpath(DATA, "fig7_weightgrid_ks.csv"), ',')
    alph_values = 2.0 .^ (-6:10)
    bet_values = 2.0 .^ (-10:7)
    p2 = heatmap(bet_values, alph_values, weightgrid, xscale=:log2, yscale=:log2, clims=(0.1, 1),
                 xlabel="Weight on Cut Edges", ylabel="Weight on Percents")

    savefig(plot(p1, p2, layout=(1, 2), size=(1000, 400)), joinpath(OUT, "figure7.png"))
    println("  wrote figures/figure7.png")
end
GC.gc()

println("=== Figure 8: Texas marginal distributions ===")
let
    cuts_rows = load_binary_gz(joinpath(DATA, "fig8_tx_cuts.bin.gz"), Int16)
    p1 = plot(xlabel="Cut Edges", ylabel="Density", xlims=(3340, 3420))
    for (i, row) in enumerate(cuts_rows)
        histogram!(p1, row, bins=3340:1:3420, normalize=true, fill=true, fillcolor=:lightblue,
                   linecolor=:blue, alpha=0.3, linealpha=0.5, label=(i == 1 ? "MEW (10 Runs)" : ""))
    end
    plot!(p1, 3340:0.5:3420, gaussian.(3340:0.5:3420, 0.01, 3346), color=:red, label="Target")
    cuts_rows = nothing; GC.gc()

    mm_rows = load_binary_gz(joinpath(DATA, "fig8_tx_mm.bin.gz"), Float32)
    p2 = plot(xlabel="Mean Median Score", ylabel="Density", xlims=(-0.010, 0.010))
    for (i, row) in enumerate(mm_rows)
        histogram!(p2, row, bins=200, normalize=true, fill=true, fillcolor=:lightblue,
                   linecolor=:blue, alpha=0.3, linealpha=0.5, label=(i == 1 ? "MEW (10 Runs)" : ""))
    end
    plot!(p2, -0.01:0.0001:0.01, gaussian.(-0.01:0.0001:0.01, 100000, 0), color=:red, label="Target")
    mm_rows = nothing; GC.gc()

    savefig(plot(p1, p2, layout=(1, 2), size=(1000, 400)), joinpath(OUT, "figure8.png"))
    println("  wrote figures/figure8.png")
end

println("Done. See $(OUT)/")

# Converts Figure 3's two partition lists into the same ragged-row-safe
# binary format as everything else in this folder (see binconv_utils.jl /
# load_data.jl): Int64 nrows, then nrows x Int64 row-lengths, then each row's
# Int8 data back-to-back. Both saved in 0/1 convention (node -> part).
#
# fig3_enumerated: all 34,225 balanced 2-partitions of Cheshire County
#   (enumpart output, space-separated 0/1 per line).
# fig3_sampled: the unique partitions actually visited by a 10M-step MEW
#   chain targeting the uniform distribution (33,750 of 34,225, 98.6%
#   coverage) -- not the 10M raw draws, just the deduplicated set, since
#   that's all the comparison-to-enumeration figure needs.

using Serialization

const NH_SRC = joinpath(@__DIR__, "..", "..", "Marked_edges")

function write_binary_gz(bin_path::String, rows::Vector{Vector{Int8}})
    open(bin_path, "w") do io
        write(io, Int64(length(rows)))
        for row in rows
            write(io, Int64(length(row)))
        end
        for row in rows
            write(io, row)
        end
    end
    run(`gzip -f -9 $bin_path`)
    gz_path = bin_path * ".gz"
    println("  wrote $(basename(gz_path)) ($(filesize(gz_path)) bytes)")
end

println("Converting enum_files/output_partitions.dat -> fig3_enumerated.bin.gz ...")
enum_rows = Vector{Vector{Int8}}()
open(joinpath(NH_SRC, "enum_files", "output_partitions.dat"), "r") do file
    for line in eachline(file)
        push!(enum_rows, parse.(Int8, split(line)))
    end
end
println("  $(length(enum_rows)) enumerated partitions")
write_binary_gz(joinpath(@__DIR__, "fig3_enumerated.bin"), enum_rows)

println("Converting cheshire_visited_partitions_10M.jls -> fig3_sampled.bin.gz ...")
visited = deserialize(joinpath(NH_SRC, "cheshire_visited_partitions_10M.jls"))
sampled_rows = [Int8.(collect(t) .- 1) for t in visited]  # 1/2 -> 0/1, matching enumerated
println("  $(length(sampled_rows)) unique sampled partitions ($(round(100*length(sampled_rows)/length(enum_rows), digits=1))% coverage)")
write_binary_gz(joinpath(@__DIR__, "fig3_sampled.bin"), sampled_rows)

# Shared by convert_to_binary.jl (NH) and convert_texas_to_binary.jl (TX).
#
# Binary layout (ragged-row-safe):
#   Int64 nrows
#   nrows x Int64   (length of each row)
#   then row 1's data (T, length[1] values), row 2's data (T, length[2] values), ...
#
# Rows are NOT assumed equal length: some MEW chains stopped early (e.g.
# Figure 8's Texas cuts.csv/mm.csv have one chain at ~2.57M steps instead of
# 4M). No padding/truncation -- every real sample collected is kept.
#
# Streams one CSV line at a time (eachline, not readdlm/readlines) and keeps
# only the small parsed-to-target-type row vectors across iterations, not the
# raw text -- this machine has 8GB RAM and some source CSVs are 1GB+ of text.

function convert_csv_to_binary(csv_path::String, bin_path::String, T::Type, expected_nrows::Int)
    println("Converting $(basename(csv_path)) -> $(basename(bin_path)) ($T, expecting $expected_nrows rows) ...")

    rows = Vector{Vector{T}}()
    for line in eachline(csv_path)
        fields = split(line, ',')
        row = Vector{T}(undef, length(fields))
        @inbounds for j in eachindex(fields)
            row[j] = parse(T, fields[j])
        end
        push!(rows, row)
        println("  row $(length(rows)): len=$(length(row)) min=$(minimum(row)) max=$(maximum(row))")
    end
    length(rows) == expected_nrows || error("$(basename(csv_path)): expected $expected_nrows rows, got $(length(rows))")

    open(bin_path, "w") do io
        write(io, Int64(length(rows)))
        for row in rows
            write(io, Int64(length(row)))
        end
        for row in rows
            write(io, row)
        end
    end
    println("  wrote $(basename(bin_path)) ($(filesize(bin_path)) bytes)")

    run(`gzip -f -9 $bin_path`)
    gz_path = bin_path * ".gz"
    println("  gzipped -> $(basename(gz_path)) ($(filesize(gz_path)) bytes)\n")
end

# Reads the .bin.gz files written by convert_to_binary.jl / convert_texas_to_binary.jl.
# Format: Int64 nrows, then nrows x Int64 row-lengths, then each row's data of
# type `T` back-to-back (row i = chain/run i). Rows may differ in length (a
# chain that stopped early) -- this always returns a Vector{Vector{T}}, never
# a Matrix, so ragged data doesn't need special-casing by callers.
# Uses the system `gunzip`, no package deps.

function load_binary_gz(path::String, ::Type{T}) where T
    open(`gunzip -c $path`, "r") do io
        nrows = read(io, Int64)
        lens = [read(io, Int64) for _ in 1:nrows]
        rows = Vector{Vector{T}}(undef, nrows)
        for i in 1:nrows
            row = Vector{T}(undef, lens[i])
            read!(io, row)
            rows[i] = row
        end
        return rows
    end
end

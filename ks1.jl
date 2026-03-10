using DelimitedFiles, Plots, ProgressBars, Distributed
addprocs(45)

@everywhere begin
    using DelimitedFiles, ProgressBars
    
    function build_count_grid_fast!(grid, data, xs, ys, x_dict, y_dict)
        fill!(grid, 0)
        for i in 1:size(data, 1)
            xi = get(x_dict, data[i,1], 0)
            yi = get(y_dict, data[i,2], 0)
            if xi > 0 && yi > 0
                grid[xi, yi] += 1
            end
        end
        
        # Cumulative sum in-place
        for i in 2:size(grid, 1)
            @views grid[i, :] .+= grid[i-1, :]
        end
        for j in 2:size(grid, 2)
            @views grid[:, j] .+= grid[:, j-1]
        end
    end
    
    function ks_two_sample_2d_fast(data1, data2, grid1, grid2, xs, ys, x_dict, y_dict)
        build_count_grid_fast!(grid1, data1, xs, ys, x_dict, y_dict)
        build_count_grid_fast!(grid2, data2, xs, ys, x_dict, y_dict)
        
        n1 = size(data1, 1)
        n2 = size(data2, 1)
        n1_inv = 1.0 / n1
        n2_inv = 1.0 / n2
        
        maxdist = 0.0
        for i in 1:length(xs)
            for j in 1:length(ys)
                d = abs(grid1[i, j] * n1_inv - grid2[i, j] * n2_inv)
                maxdist = max(maxdist, d)
            end
        end
        return maxdist
    end
    
    function incremental_ks_2samp_2d_fast(data1, data2)
        max_iterations = min(size(data1, 1), size(data2, 1))
        
        # Pre-allocate everything based on full dataset to avoid repeated allocations
        all_data = vcat(data1, data2)
        xs = sort!(unique(all_data[:,1]))
        ys = sort!(unique(all_data[:,2]))
        
        # Create lookup dictionaries for O(1) index finding
        x_dict = Dict(v => i for (i, v) in enumerate(xs))
        y_dict = Dict(v => i for (i, v) in enumerate(ys))
        
        # Pre-allocate grids
        grid1 = zeros(Int, length(xs), length(ys))
        grid2 = zeros(Int, length(xs), length(ys))
        
        # Pre-allocate result arrays
        kss = Float64[]
        iterations = Int[]
        
        # For consistent loglog spacing
        next_sample = 1
        multiplier = 1.9
        
        for i in ProgressBar(1:max_iterations)
            if i >= next_sample
                # Calculate 2D KS distance using views (no copying)
                current_data1 = @view data1[1:i, :]
                current_data2 = @view data2[1:i, :]
                
                ks_dist = ks_two_sample_2d_fast(current_data1, current_data2, grid1, grid2, xs, ys, x_dict, y_dict)
                push!(kss, ks_dist)
                push!(iterations, i)
                next_sample = max(i + 1, round(Int, next_sample * multiplier))
            end
        end
        return iterations, kss
    end
    
    function wrapper_precomputed(params, data_pairs)
        data1, data2 = data_pairs[(params.i, params.j)]
        return incremental_ks_2samp_2d_fast(data1, data2)
    end
end

# Pre-compute only the data combinations we need
function prepare_data_pairs()
    println("Loading data files...")
    cuts = readdlm("cuts.csv", ',')
    mms = readdlm("mm.csv", ',')
    
    println("Preparing data pairs...")
    data_pairs = Dict()
    for i in 1:10, j in i+1:10
        key = (i, j)
        data1 = Matrix{Float64}(hcat(cuts[i,:], mms[i,:]))
        data2 = Matrix{Float64}(hcat(cuts[j,:], mms[j,:]))
        data_pairs[key] = (data1, data2)
    end
    println("Data preparation complete. Created $(length(data_pairs)) pairs.")
    return data_pairs
end

# Main execution
println("Starting computation...")
data_pairs = prepare_data_pairs()
params = [(i=i, j=j) for i in 1:10 for j in i+1:10]

println("Processing $(length(params)) parameter combinations...")
ksss = pmap(p -> wrapper_precomputed(p, data_pairs), params)

println("Writing results...")
writedlm("ksss_2d.csv", ksss, ',')
println("Complete! Results saved to ksss_2d.csv")
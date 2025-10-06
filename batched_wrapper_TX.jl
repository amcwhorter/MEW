using Distributed

addprocs(10) 

@everywhere begin
    include("main_TX.jl")s

    function long_run(n, j, alph, bet, m1, m2)
        GC.gc()
        run_dir = "push_runs1/run$(j)"
        if isdir(run_dir)
            existing_files = filter(f -> endswith(f, ".jls"), readdir(run_dir))
            if !isempty(existing_files)
                run_numbers = [parse(Int, match(r"run(\d+)\.jls", f).captures[1]) for f in existing_files]
                last_run = maximum(run_numbers)
                println("Found existing runs up to run$(last_run) in $(run_dir)")
                last_file = "$(run_dir)/run$(last_run).jls"
                c_last = deserialize(last_file)
                initialization = c_last[2:3]
                start_i = last_run + 1
            else
                start_i = 1
                c1 = main(alph, bet, m1, m2)
                initialization = c1[2:3]
            end
        else
            println("Creating new directory: $(run_dir)")
            mkpath(run_dir)
            c1 = main(alph, bet, m1, m2)
            initialization = c1[2:3]
            start_i = 1
        end
        if start_i <= n
            for i in start_i:n
                c2 = main(alph, bet, m1, m2, initialization)
                serialize("$(run_dir)/run$(i).jls", c2)
                println("starting run $i on run $j")
                initialization = c2[2:3]
            end
        else
            println("All $(n) runs already completed for this parameter set.")
        end
        return nothing

    end

    function run_wrapper(params)
        j = params.j
        bet = params.bet
        m2 = params.m2
        println("Starting run $j with m2 = $m2, beta=$bet")
        # try 
        long_run(num_iterations, j, alph, bet, m1, m2)
        # catch e
        #     println("Error in run $j with alpha = $alph, beta = $bet: $e")
        # end
        return nothing
    end

    num_iterations = 10000
    m1 = 1
    m2_values = [2]

    alph = 10
    bet_values = [1]



    j_values = 1:10

    param_combinations = [(j=j, bet=bet, m2=m2) for j in j_values, bet in bet_values, m2 in m2_values]
    param_combinations = vec(param_combinations)
end



# run_wrapper(param_combinations[1])

results = pmap(run_wrapper, param_combinations)

using Dierckx

# Structure Refinement
function Refine(data::Matrix{T}, tmatch::Matrix{T}, ensize::Int64, qstart::T; population::Int64=10, maxiter::Int64=100000, reports::Int64=100) where T<:Real

    # Initialize
    local NGENE = size(tmatch, 2) - 1;
    local GENES = collect(1:NGENE);
    chromosomes = reshape(sample(GENES, ensize*population, replace=true), ensize, population);


    # Pre-process the data according to qstart
    id = findall(data[:,1] .> qstart)[1];
    expdata = data[id:end, :];

    # filter out Problematic data
    expdata = expdata[(expdata[:, 3] .!= 0), :]; # error = 0.0
    expdata[:, 2] = abs.(expdata[:, 2]);

    t = zeros(size(expdata,1), size(tmatch,2));
    [t[:,i] = Spline1D(tmatch[:, 1], tmatch[:, i]; k=1)(expdata[:,1]) for i = 1:size(tmatch, 2)];
    t[:,1] = expdata[:,1];


    # Fitness Value
    function fitness(data::Matrix{T}, t::Matrix{T}, chromosomes::Matrix{Int64}, population::Int64) where T<:Real
        ft = zeros(population);

        function logfit(iexp::Matrix{T}, ifit::Vector{T}) where T<:Real

            #σh = iexp[:,2] + iexp[:,3];
            #σl = iexp[:,2] - iexp[:,3];

            σ = iexp[:,3]./(iexp[:,2] * log(10));
            #σ[index_pos] = 0.5 * (log10.(σh[index_pos]) - log10.(σl[index_pos]));
            return mean(((log10.(iexp[:,2]) - log10.(ifit))./σ).^2);
        end

        [ft[i] = logfit(expdata, mean(t[:, chromosomes[:,i] .+ 1], dims=2)[:,1]) for i = 1:population];
        return ft;
    end


    # Selection
    function selection(chromosomes::Matrix{Int64}, ft::Vector{T}, population::Int64) where T<:Real
        m = Int64(floor(population/3));
        sele = sample(collect(1:population), ProbabilityWeights(exp.(-3.0*ft)), m, replace=false);
        parents = chromosomes[:, sele];
        return parents;
    end

    # Crossover
    function crossover(parents::Matrix{Int64}, ensize::Int64, population::Int64)
        child = Matrix{Int64}(undef, ensize, population);
        nparents = size(parents, 2);
        for i = 1:population
            p = sample(collect(1:nparents), 2);
            father, mother = parents[:, p[1]], parents[:, p[2]];
            sameindex = findall(indexin(father, mother) .!== nothing);
            m = length(sameindex);
            if m == 0
                child[:, i] = sample([father; mother], ensize, replace=false);
            else
                child[1:m, i] = father[sameindex];
                child[m+1:end, i] = sample(unique([father; mother]), ensize-m);
            end
        end
        return child;
    end

    # Mutation
    function mutation(child::Matrix{Int64}, ensize::Int64, population::Int64; rate::Float64=0.3)
        mutated_child = copy(child);
        for i = 1:population
            if rand() <= rate
                nmutation = sample(collect(1:ensize), ProbabilityWeights(1.0 ./collect(1:ensize)));
                mutated_child[sample(collect(1:ensize), nmutation, replace=false), i] = sample(GENES, nmutation, replace=false);
            end
        end
        return mutated_child;
    end

    # Main process
    best_fitness = Inf;
    best_chromosome = Vector{Int64}(undef, size(chromosomes, 1));

    tick = Int64(round(maxiter/reports));

    for i = 1:maxiter
        ft = fitness(expdata, t, chromosomes, population);
        parents = selection(chromosomes, ft, population);
        child = crossover(parents, ensize, population);
        mutated_child = mutation(child, ensize, population);

        m = minimum(ft);

        if m <= best_fitness
            best_fitness = m;
            best_chromosome = chromosomes[:, findall(ft.==m)];
        end

        chromosomes = mutated_child;
        i % tick == 0 ? println("Generation #$i: fitness = $(round(m, digits=5)), best fit so far = $(round(best_fitness, digits=5)).") : nothing;
    end
    println("Finished Genetic Algorithm, with best fitness = $(round(best_fitness, digits=5)).");
    return best_fitness, best_chromosome;
end

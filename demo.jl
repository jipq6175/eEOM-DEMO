
# Julia script for eEOM demo
# Author: Yen-Lin Chen


# Load the dependencies
using Plots, DelimitedFiles, StatsBase, JLD
include("ga.jl");

# Load data
pool1, pool2 = load("pool.jld", "pool1", "pool2");
data1, data2 = readdlm("data1.dat"), readdlm("data2.dat");
q = collect(0.0:0.005:1.25);


## First look at the [data1, pool1] set

# Plot the data and see what's going on
p = plot(size=(400, 300), dpi=600, grid=false, legend=:outerbottomright, framestyle=:box);
s = 0.4;
for i = 1:14
    plot!(data1[:,1], s^(i-1) * data1[:,2], yscale=:log10, size=(800,900), lw=2.0, c=:black, linealpha=0.5, lab="");
    plot!(q, s^(i-1)*pool1[:,7i-6], lw=1.5, yscale=:log10, lab="Conf. $(12+i)");
end
xlims!(0.05, 0.8)
ylims!(0.8e-6, 10)
ylabel!("I(q) (A.U.)")
xlabel!("q (1/A)")


# Run the eEOM structural refinement
idx = findall(data1[:,1] .<= 0.8);
tmp = data1[idx, :];
ensize = 10;
cycles = 50;
survivors1 = zeros(Int64, ensize, cycles);
fitness1 = zeros(cycles);

t = @elapsed for i = 1:cycles
    @info("===== Cycle number $i ====");
    fitness1[i], survivor = Refine(tmp, [q pool1], ensize, 0.15; population=10, maxiter=100000, reports=10);
    survivors1[:,i] = survivor[:,1];
end

println("- eEOM for $cycles cycles and $(size(pool1, 2)) conformaitons took $(round(t, digits=2)) seconds.");



# Plot the eEOM results
p = scatter(data1[:,1], data1[:,2], yerror=data1[:,3], yscale=:log10, markersize=5.0, c=:black, markeralpha=0.6, linealpha=0.25, lab="Experiment", size=(300,450), legendfont=font(16), titlefont=font(20));
plot!(q, mean(pool1[:, survivors1[:, findall(fitness1 .== minimum(fitness1))[1]]], dims=2), lw=5, linealpha=0.8, lab="eEOM Fit")
plot!(size=(600, 800), dpi=600, grid=false, legend=true, framestyle=:box);
xlims!(0.05, 0.8);
ylims!(0.1, 10)
xlabel!("q (1/A)");
ylabel!("I(q) (A.U.)")


# Plot the best ensemble
best = survivors1[:, findall(fitness1 .== minimum(fitness1))[1]];
dic = countmap(reshape(best, (length(best), )));
b = bar(collect(keys(dic)), collect(values(dic)) / length(best), title="Cluster Selection: Best", xticks=collect(1:2:98), bar_width=0.5);
plot!(size=(500, 300), dpi=600, grid=false, legend=false, framestyle=:box);
ylabel!("Population Fraction")


# Plot the all ensemble
dic = countmap(reshape(survivors1, (length(survivors1), )));
b = bar(collect(keys(dic)), collect(values(dic)) / length(best), title="Cluster Selection: All", xticks=collect(1:2:98), bar_width=0.5);
plot!(size=(500, 300), dpi=600, grid=false, legend=false, framestyle=:box);
ylabel!("Population Fraction")


# And then do whatever you like from the results


## Second look at the [data2, pool2] set

# Plot the data and see what's going on
p = plot(size=(400, 300), dpi=600, grid=false, legend=:outerbottomright, framestyle=:box);
s = 0.4;
for i = 1:14
    plot!(data2[:,1], s^(i-1) * data2[:,2], yscale=:log10, size=(800,900), lw=2.0, c=:black, linealpha=0.5, lab="");
    plot!(q, s^(i-1)*pool2[:,7i-6], lw=1.5, yscale=:log10, lab="Conf. $(12+i)");
end
xlims!(0.05, 0.8)
ylims!(0.8e-6, 10)
ylabel!("I(q) (A.U.)")
xlabel!("q (1/A)")


# Run the eEOM structural refinement
idx = findall(data2[:,1] .<= 1.0);
tmp = data2[idx, :];
ensize = 10;
cycles = 50;
survivors2 = zeros(Int64, ensize, cycles);
fitness2 = zeros(cycles);

t = @elapsed for i = 1:cycles
    @info("===== Cycle number $i ====");
    fitness2[i], survivor = Refine(tmp, [q pool2], ensize, 0.10; population=10, maxiter=100000, reports=10);
    survivors2[:,i] = survivor[:,1];
end

println("- eEOM for $cycles cycles and $(size(pool2, 2)) conformaitons took $(round(t, digits=2)) seconds.");



# Plot the eEOM results
p = scatter(data2[:,1], data2[:,2], yerror=data2[:,3], yscale=:log10, markersize=5.0, c=:black, markeralpha=0.6, linealpha=0.25, lab="Experiment", size=(300,450), legendfont=font(16), titlefont=font(20));
plot!(q, mean(pool2[:, survivors2[:, findall(fitness2 .== minimum(fitness2))[1]]], dims=2), lw=5, linealpha=0.8, lab="eEOM Fit")
plot!(size=(600, 800), dpi=600, grid=false, legend=true, framestyle=:box);
xlims!(0.05, 0.8);
ylims!(0.1, 10)
xlabel!("q (1/A)");
ylabel!("I(q) (A.U.)")


# Plot the best ensemble
best = survivors2[:, findall(fitness2 .== minimum(fitness2))[1]];
dic = countmap(reshape(best, (length(best), )));
b = bar(collect(keys(dic)), collect(values(dic)) / length(best), title="Cluster Selection: Best", xticks=collect(1:4:98), bar_width=0.5);
plot!(size=(500, 300), dpi=600, grid=false, legend=false, framestyle=:box);
ylabel!("Population Fraction")


# Plot the all ensemble
dic = countmap(reshape(survivors2, (length(survivors2), )));
b = bar(collect(keys(dic)), collect(values(dic)) / length(best), title="Cluster Selection: All", xticks=collect(1:4:98), bar_width=0.5);
plot!(size=(500, 300), dpi=600, grid=false, legend=false, framestyle=:box);
ylabel!("Population Fraction")


# And then do whatever you like from the results


# EOF

using Distributions

n_sys = 4
n_comp = 5
# set parameters
lams = zeros(n_comp, n_sys)

for i in 1:n_comp
        lams[i,:] = .1 * i * ones(n_sys)
end

testLen = 500

df = DataFrame(A = Int64[], B = Int64[])

julia> push!(df, [3  6])

for sys in 1:n_sys
        for comp in 1:n_comp
                t1 = 0

                while t1 < testLen

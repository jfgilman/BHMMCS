using Distributions
using DataFrames
using Mamba

n_sys = 4
n_comp = 5
# set parameters
lams = zeros(n_comp, n_sys)

for i in 1:n_comp
        lams[i,:] = .1 * i * ones(n_sys)
end

testLen = 500

df = DataFrame(Time = Float64[], Censored = Bool[], System = Int64[])

push!(df, [3.14, true, 1])

for sys in 1:n_sys
        for comp in 1:n_comp
                t1 = 0

                while t1 < testLen

                end
        end
end


## Data
JLTV = Dict{Symbol, Any}(
  :phase => [1, 1, 1, 2, 2, 2, 3, 3, 3],
  :system => [1, 1, 1, 1, 2, 2, 2, 2, 2],
  :mbf => [1.3, 3.2, 7.3, 3.1, 5.2, 3.4, 2.0, 8.8, 10.2],
)
JLTV[:n_sys] = 4
JLTV[:N] = length(JLTV[:mbf])

model = Model(

  mbf = Stochastic(1,
      (lambda, phase, system, N, rho1, rho2) ->
        UnivariateDistribution[
          begin
            if phase[i] == 1
              rate = lambda[system[i]]
            elseif phase[i] == 2
              rate = lambda[system[i]] * rho1
            else
              rate = lambda[system[i]] * rho1 * rho2
            end
            Exponential(rate)
          end
          for i in 1:N
        ],
      false
  ),

  lambda = Stochastic(1,
    (alpha, beta) -> Gamma(alpha, 1 / beta),
    true
  ),

  alpha = Stochastic(
    () -> Gamma(0.001, 0.001)
  ),

  beta = Stochastic(
    () -> Gamma(0.001, 0.001)
  ),

  rho1 = Stochastic(
    () -> Gamma(3.0, 3.0)
  ),

  rho2 = Stochastic(
    () -> Gamma(3.0, 3.0)
  )
)

## Initial Values
inits = [
  Dict{Symbol, Any}(
    :mbf => JLTV[:mbf],
    :lambda => rand(Gamma(1, 1), JLTV[:n_sys]),
    :alpha => rand(Gamma(1, 1)),
    :beta => rand(Gamma(1, 1)),
    :rho1 => rand(Gamma(1, 1)),
    :rho2 => rand(Gamma(1, 1))
  )
for i in 1:3
]

## No-U-Turn Sampling Scheme
scheme2 = [NUTS([:lambda, :alpha, :beta, :rho1, :rho2])]

## Sampling Scheme Assignment
setsamplers!(model, scheme2)
sim2 = mcmc(model, JLTV, inits, 10000, burnin=250, thin=2, chains=3)

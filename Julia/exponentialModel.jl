using Distributions
using DataFrames
using Mamba

n_sys = 4
n_comp = 5
# set parameters
lams = zeros(n_comp, n_sys)

for i in 1:n_sys
        lams[:,i] = 5 * i * ones(n_comp)
end

rho1 = 0.7
rho2 = 0.9

testLen = 1000

df = DataFrame(MBF = Float64[], Censored = Bool[], System = Int64[], Phase = Int64[])

for sys in 1:n_sys
      t1 = 0
      while t1 < testLen
        d = rand(Exponential(lams[1, sys]))
        t1 += d

        if t1 > testLen
          push!(df, [d, true, sys, 1])
        else
          push!(df, [d, false, sys, 1])
        end
      end

      t2 = 0
      while t2 < testLen
        d = rand(Exponential(lams[1, sys] * rho1))
        t2 += d

        if t2 > testLen
          push!(df, [d, true, sys, 2])
        else
          push!(df, [d, false, sys, 2])
        end
      end

      t3 = 0
      while t3 < testLen
        d = rand(Exponential(lams[1, sys] * rho1 * rho2))
        t3 += d

        if t3 > testLen
          push!(df, [d, true, sys, 3])
        else
          push!(df, [d, false, sys, 3])
        end
      end
end

## Data
JLTV = Dict{Symbol, Any}(
  :phase => Array(df[:Phase]),
  :system => Array(df[:System]),
  :mbf => Array(df[:MBF]),
  :censored => Array(df[:Censored])
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
    (alpha, beta) -> Gamma(alpha, beta),
    true
  ),

  alpha = Stochastic(
    () -> Gamma(0.001, 0.001),
    false
  ),

  beta = Stochastic(
    () -> Gamma(0.001, 0.001),
    false
  ),

  rho1 = Stochastic(
    () -> Gamma(3.0, 3.0)
  ),

  rho2 = Stochastic(
    () -> Gamma(3.0, 3.0)
  )
)

n_chains = 2
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
for i in 1:n_chains
]

Gibbs_beta = Sampler([:beta],
  (beta, alpha, lambda) ->
    begin
      hyper_prior = 0.001
      a = 4 + shape(beta.distr)
      sum_lam = sum(lambda)
      b = sum_lam + scale(beta.distr)
      rand(Gamma(a, b))
    end
)

# Gibbs_lambda = Sampler([:lambda],
#   (beta, alpha, lambda, rho1, rho2, mbf, phase, censored) ->
#     begin
#       a = alpha
#
#       rand(Gamma(a, b))
#     end
# )

scheme1 = [Gibbs_beta, NUTS([:lambda, :alpha, :rho1, :rho2])]

## No-U-Turn Sampling Scheme
scheme2 = [NUTS([:lambda, :alpha, :beta, :rho1, :rho2])]

## Sampling Scheme Assignment
setsamplers!(model, scheme1)
sim2 = mcmc(model, JLTV, inits, 10000, burnin=250, thin=2, chains=n_chains)

describe(sim2)

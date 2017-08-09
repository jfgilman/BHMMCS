using Distributions
using DataFrames
using Mamba

n_sys = 4
n_comp = 5
# set parameters
lams = zeros(n_comp, n_sys)

for i in 1:n_comp
        lams[i,:] = 5 * i * ones(n_sys)
end

rho1 = 0.7
rho2 = 0.9

testLen = 500

df = DataFrame(MBF = Float64[], Censored = Bool[], System = Int64[], Phase = Int64[])

for sys in 1:n_sys
        for comp in 1:n_comp
                t1 = 0

                while t1 < testLen
                  d = rand(Exponential(lams[comp, sys]))
                  t1 += d

                  if t1 > testLen
                    push!(df, [d, true, sys, 1])
                  else
                    push!(df, [d, false, sys, 1])
                  end
                end

                t2 = 0
                while t2 < testLen
                  d = rand(Exponential(lams[comp, sys] * rho1))
                  t2 += d

                  if t2 > testLen
                    push!(df, [d, true, sys, 2])
                  else
                    push!(df, [d, false, sys, 2])
                  end
                end

                t3 = 0
                while t3 < testLen
                  d = rand(Exponential(lams[comp, sys] * rho1 * rho2))
                  t3 += d

                  if t3 > testLen
                    push!(df, [d, true, sys, 3])
                  else
                    push!(df, [d, false, sys, 3])
                  end
                end
        end
end

## Data
JLTV = Dict{Symbol, Any}(
  :phase => Array(df[:Phase]),
  :system => Array(df[:System]),
  :mbf => Array(df[:MBF])
)
JLTV[:n_sys] = 4
JLTV[:N] = length(JLTV[:mbf])

# ## Data
# JLTV = Dict{Symbol, Any}(
#   :phase => [1, 1, 2, 2, 3, 3, 1],
#   :system => [1,2,3,4,1,2,3],
#   :mbf => [2.0, 3.3332, 2.22, 5.4421, 0.2934, 10.2, 11.0]
# )
# JLTV[:n_sys] = 4
# JLTV[:N] = length(JLTV[:mbf])

typeof(JLTV[:mbf])

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
            Exponential(1 / rate)
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

Gibbs_beta = Sampler([:beta],
  (beta, s2, xmat, y) ->
    begin
      beta_mean = mean(beta.distr)
      beta_invcov = invcov(beta.distr)
      Sigma = inv(Symmetric(xmat' * xmat / s2 + beta_invcov))
      mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
      rand(MvNormal(mu, Sigma))
    end
)

## No-U-Turn Sampling Scheme
scheme2 = [NUTS([:lambda, :alpha, :beta, :rho1, :rho2])]

## Sampling Scheme Assignment
setsamplers!(model, scheme2)
sim2 = mcmc(model, JLTV, inits, 10000, burnin=250, thin=2, chains=3)


describe(sim2)

using Mamba

## Model Specification

model = Model(

  y = Stochastic(1,
    (mu, s2) ->  MvNormal(mu, sqrt(s2)),
    false
  ),

  mu = Logical(1,
    (xmat, beta) -> xmat * beta,
    false
  ),

  beta = Stochastic(1,
    () -> MvNormal(2, sqrt(1000))
  ),

  s2 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  )

)

# # Case 1: Multivariate Normal with independence covariance matrix
# beta = Stochastic(1,
#   () -> MvNormal(2, sqrt(1000))
# )
#
# # Case 2: One common univariate Normal
# beta = Stochastic(1,
#   () -> Normal(0, sqrt(1000))
# )
#
# # Case 3: Array of univariate Normals
# beta = Stochastic(1,
#   () -> UnivariateDistribution[Normal(0, sqrt(1000)), Normal(0, sqrt(1000))]
# )
#
# # Case 4: Array of univariate Normals
# beta = Stochastic(1,
#   () -> UnivariateDistribution[Normal(0, sqrt(1000)) for i in 1:2]
# )
#
# ## Hybrid No-U-Turn and Slice Sampling Scheme
# scheme1 = [NUTS(:beta),
#            Slice(:s2, 3.0)]
#
# ## No-U-Turn Sampling Scheme
# scheme2 = [NUTS([:beta, :s2])]


## User-Defined Samplers

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

Gibbs_s2 = Sampler([:s2],
  (mu, s2, y) ->
    begin
      a = length(y) / 2.0 + shape(s2.distr)
      b = sumabs2(y - mu) / 2.0 + scale(s2.distr)
      rand(InverseGamma(a, b))
    end
)

## User-Defined Sampling Scheme
scheme3 = [Gibbs_beta, Gibbs_s2]

## Data
line = Dict{Symbol, Any}(
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
)
line[:xmat] = [ones(5) line[:x]]

## Initial Values
inits = [
  Dict{Symbol, Any}(
    :y => line[:y],
    :beta => rand(Normal(0, 1), 2),
    :s2 => rand(Gamma(1, 1))
  )
for i in 1:3
]

setsamplers!(model, scheme3)
sim3 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)

gelmandiag(sim3, mpsrf=true, transform=true) |> showall

heideldiag(sim3) |> showall

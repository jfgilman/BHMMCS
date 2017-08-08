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


################################################################################
# Example 2
################################################################################



## Data
birats = Dict{Symbol, Any}(
  :N => 30, :T => 5,
  :x => [8.0, 15.0, 22.0, 29.0, 36.0],
  :Y => [151 199 246 283 320
         145 199 249 293 354
         147 214 263 312 328
         155 200 237 272 297
         135 188 230 280 323
         159 210 252 298 331
         141 189 231 275 305
         159 201 248 297 338
         177 236 285 350 376
         134 182 220 260 296
         160 208 261 313 352
         143 188 220 273 314
         154 200 244 289 325
         171 221 270 326 358
         163 216 242 281 312
         160 207 248 288 324
         142 187 234 280 316
         156 203 243 283 317
         157 212 259 307 336
         152 203 246 286 321
         154 205 253 298 334
         139 190 225 267 302
         146 191 229 272 302
         157 211 250 285 323
         132 185 237 286 331
         160 207 257 303 345
         169 216 261 295 333
         157 205 248 289 316
         137 180 219 258 291
         153 200 244 286 324],
  :mean => [0.0, 0.0],
  :var => [1.0e6 0.0
           0.0   1.0e6],
  :Omega => [200.0 0.0
               0.0 0.2]
)


## Model Specification
model = Model(

  Y = Stochastic(2,
    (beta, x, sigmaC, N, T) ->
      UnivariateDistribution[
        Normal(beta[i, 1] + beta[i, 2] * x[j], sigmaC)
        for i in 1:N, j in 1:T
      ],
    false
  ),

  beta = Stochastic(2,
    (mu_beta, Sigma, N) ->
      MultivariateDistribution[
        MvNormal(mu_beta, Sigma)
        for i in 1:N
      ],
    false
  ),

  mu_beta = Stochastic(1,
    (mean, var) -> MvNormal(mean, var)
  ),

  Sigma = Stochastic(2,
    Omega -> InverseWishart(2, Omega),
    false
  ),

  sigma2C = Stochastic(
    () -> InverseGamma(0.001, 0.001),
    false
  ),

  sigmaC = Logical(
    sigma2C -> sqrt(sigma2C)
  )

)


## Initial Values
inits = [
  Dict(:Y => birats[:Y], :beta => repmat([100 6], birats[:N], 1),
       :mu_beta => [0, 0], :Sigma => eye(2), :sigma2C => 1.0),
  Dict(:Y => birats[:Y], :beta => repmat([50 3], birats[:N], 1),
       :mu_beta => [10, 10], :Sigma => 0.3 * eye(2), :sigma2C => 10.0)
]


## Sampling Scheme
scheme = [AMWG([:beta, :mu_beta], repmat([10.0, 1.0], birats[:N] + 1)),
          AMWG(:Sigma, 1.0),
          Slice(:sigma2C, 10.0)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, birats, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)

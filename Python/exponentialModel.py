"""Exponential Model."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc3 as pm
import theano.tensor as t
from theano import function

plt.style.use('ggplot')

df = pd.read_csv('Python/simDat.csv')

df1 = df[df.System == 1]

n_component = max(df['Component']).astype(int)

times = df1.Time
fail = 1 - (df1.Censored).astype(int)
comp = (df1.Component - 1).astype(int)
p2 = df1.Phase2.values
p3 = df1.Phase3.values

lam1 = (df.Component[:, None] == 1).astype(int)
lam2 = (df.Component[:, None] == 2).astype(int)
lam3 = (df.Component[:, None] == 3).astype(int)
lam4 = (df.Component[:, None] == 4).astype(int)
lam5 = (df.Component[:, None] == 5).astype(int)

mini_batch_size = 100

v_lamIn1 = pm.Minibatch(lam1, batch_size=mini_batch_size)
v_lamIn2 = pm.Minibatch(lam2, batch_size=mini_batch_size)
v_lamIn3 = pm.Minibatch(lam3, batch_size=mini_batch_size)
v_lamIn4 = pm.Minibatch(lam4, batch_size=mini_batch_size)
v_lamIn5 = pm.Minibatch(lam5, batch_size=mini_batch_size)
v_time = pm.Minibatch(df1.Time, batch_size=mini_batch_size)
v_fail = pm.Minibatch(fail, batch_size=mini_batch_size)
v_comp = pm.Minibatch(comp, batch_size=mini_batch_size)
v_P2 = pm.Minibatch(df1.Phase2.astype(int), batch_size=mini_batch_size)
v_P3 = pm.Minibatch(df1.Phase3.astype(int), batch_size=mini_batch_size)

with pm.Model() as hm:
    # Hyperpriors
    a = pm.Gamma('alpha', alpha=0.001, beta=0.001)
    b = pm.Gamma('beta', alpha=0.001, beta=0.001)

    # parts of rate
    lam = pm.Gamma('rate', alpha=a, beta=b, shape=n_component)
    rho1 = pm.Gamma('Rho1', alpha=3, beta=3)
    rho2 = pm.Gamma('Rho2', alpha=3, beta=3)

    # rateP = [lam[int(df1.Component[i] - 1)] * rho1**df1.Phase2[i] * rho2**df1.Phase3[i] for i in range(len(df1))]
    # rateP = lam[comp] * rho1**p2 * rho2**p3
    # rateP = lam[v_comp] * rho1**v_P2 * rho2**v_P3
    # logRate = np.log(lam[comp]) + np.log(rho1) * p2 + np.log(rho2) * p3
    rateP = (lam[0] * v_lamIn1) * rho1**v_P2 * rho2**v_P3 \
        + (lam[1] * v_lamIn2) * rho1**v_P2 * rho2**v_P3 \
        + (lam[2] * v_lamIn3) * rho1**v_P2 * rho2**v_P3 \
        + (lam[3] * v_lamIn4) * rho1**v_P2 * rho2**v_P3 \
        + (lam[4] * v_lamIn5) * rho1**v_P2 * rho2**v_P3

    # rateP = (lam[0] * lam1) * rho1**p2 * rho2**p3 \
    #     + (lam[1] * lam2) * rho1**p2 * rho2**p3 \
    #     + (lam[2] * lam3) * rho1**p2 * rho2**p3 \
    #     + (lam[3] * lam4) * rho1**p2 * rho2**p3 \
    #     + (lam[4] * lam5) * rho1**p2 * rho2**p3

    def expo_log_like(failure, time):
        return t.sum(failure * t.log(rateP) - rateP * time)

    # def expo_log_like(failure, time):
    #     return t.sum(failure * logRate - t.exp(logRate) * time)

    # Data likelihood
    # y = pm.Exponential('likelihood', rateP, observed=v_time, total_size=len(times))
    y = pm.DensityDist('y', expo_log_like, observed={'failure': v_fail, 'time': v_time}, total_size=len(times))
    # y = pm.DensityDist('y', expo_log_like, observed={'failure': fail, 'time': times})

    trace = pm.sample(1000, tune=500)

pm.traceplot(trace)
plt.show()

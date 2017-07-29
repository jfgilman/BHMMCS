"""Exponential Model"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.stats as stats


plt.style.use('ggplot')

df = pd.read_csv('Python/simDat.csv')

df1 = df[df.System == 1]

n_component = max(df['Component'])

with pm.Model() as hm:
    # Hyperpriors
    a = pm.Gamma('alpha', alpha=0.001, beta=0.001)
    b = pm.Gamma('beta', alpha=0.001, beta=0.001)

    # parts of rate
    lam = pm.Gamma('rate', alpha=a, beta=b, shape=n_component)
    rho1 = pm.Gamma('Rho1', alpha=3, beta=3)
    rho2 = pm.Gamma('Rho2', alpha=3, beta=3)

    rateP = [lam[int(df1.Component[i] - 1)] * rho1**df1.Phase2[i] * rho2**df1.Phase3[i] for i in range(len(df1))]

    # Data likelihood
    y = pm.Exponential('y', rateP, observed=df1.Time)

    trace = pm.sample(1000, tune=500)

pm.traceplot(trace)
plt.show()

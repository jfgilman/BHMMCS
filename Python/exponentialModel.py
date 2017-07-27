"""Exponential Model"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.stats as stats

plt.style.use('ggplot')

df = pd.read_csv('simDat.csv')

miles = df.iloc[:20, 1]

with pm.Model() as hm:
    # County hyperpriors
    a = pm.Gamma('alpha', alpha=0.001, beta=0.001)
    b = pm.Gamma('beta', alpha=0.001, beta=0.001)

    # County slopes and intercepts
    lam = pm.Gamma('rate', alpha=a, beta=b)

    # Data likelihood
    y = pm.Exponential('y', lam=lam, observed=miles)

    trace = pm.sample(1000, tune=500)

pm.traceplot(trace)

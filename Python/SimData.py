"""Setting up sample data for 4 systems with 5 components."""

import pandas as pd
import numpy as np

# set parameters
lams = [[50] * 4,
        np.random.gamma(10, 15, 4),
        [24] * 4,
        np.random.gamma(5, 8, 4),
        [80] * 4]

rho1 = .7
rho2 = .8

testLen = 500

df = pd.DataFrame(columns=['Time', 'Censored', 'System', 'Phase', 'Component'])

i = 0
for sys in range(4):
    for comp in range(5):
        t1 = 0

        while t1 < testLen:
            d = np.random.exponential(lams[comp][sys])
            t1 += d
            if(t1 > testLen):
                df.loc[i] = [d - (t1 - testLen), True, sys + 1, 1, comp + 1]
                i += 1
            else:
                df.loc[i] = [d, False, sys + 1, 1, comp + 1]
                i += 1

        t2 = 0
        while t2 < testLen:
            d = np.random.exponential(lams[comp][sys] * rho1)
            t2 += d
            if(t2 > testLen):
                df.loc[i] = [d - (t2 - testLen), True, sys + 1, 2, comp + 1]
                i += 1
            else:
                df.loc[i] = [d, False, sys + 1, 2, comp + 1]
                i += 1

        t3 = 0
        while t3 < testLen:
            d = np.random.exponential(lams[comp][sys] * rho1 * rho2)
            t3 += d
            if(t3 > testLen):
                df.loc[i] = [d - (t3 - testLen), True, sys + 1, 3, comp + 1]
                i += 1
            else:
                df.loc[i] = [d, False, sys + 1, 3, comp + 1]
                i += 1

df.to_csv('simDat.csv')

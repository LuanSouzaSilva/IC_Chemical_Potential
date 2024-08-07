import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv(r'IC_Chemical_Potential/DMRG/DMRG_CSVs\N20U10.csv')

onsite = np.array(df['Onsite_Energy'])
mu = np.array(df['Chemical_Potential'])

onsite_ren = 1/2 + onsite/10

x = np.copy(onsite_ren)

y = x

plt.plot(x, y, color = 'black', label = r'f(x) = x')
plt.scatter(onsite_ren, mu/10, label = 'U=10', s = 50)

plt.xlabel(r'$1/2-\dfrac{\epsilon_d}{U}$', fontsize = 20)
plt.ylabel(r'$\dfrac{\mu}{U}$', fontsize = 20)

plt.axhline(0, c = 'g', linestyle = 'dashed', label = r'$\mu = 0$')

plt.tick_params(labelsize = 15)

plt.legend()

plt.tight_layout

plt.show()
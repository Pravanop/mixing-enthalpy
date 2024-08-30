import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib

matplotlib.rcParams['font.size'] = 16

df = pd.read_csv("Pravan_bcc_4_calphad.csv")
df["error"] = abs(df["Calphad"] - df["Model"])
df = df.sort_values(by = "error", ascending = False)

query = df[(df['error'].isnull())]
query = query[(query['Calphad'].isnull()) | (query['Model'].isnull())]
print(query.info())
fig, ax = plt.subplots(2, 1, figsize=(8, 12))
ax[0].scatter(df["Calphad"], df["Model"], label="Misc_T")
ax[0].plot(np.linspace(200, 3201, 50), np.linspace(200, 3201, 50),)
ax[0].set_xlabel("Calphad")
ax[0].set_ylabel("Model")
ax[0].set_xlim([200, 3200])
ax[0].set_ylim([200, 3200])

plt.subplots_adjust(hspace=0.2)
error = abs(df["Model"] - df["Calphad"])
sns.histplot(error, ax=ax[1], label="Error", kde=True)
ax[1].set_xlabel("Error in T (K)")
plt.savefig("./error_plot.png", dpi = 300)
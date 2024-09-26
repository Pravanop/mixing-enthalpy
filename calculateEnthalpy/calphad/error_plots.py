import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib

matplotlib.rcParams['font.size'] = 16

df = pd.read_csv("/Users/pravanomprakash/Documents/Projects/mixing-enthalpy/calculateEnthalpy/calphad/calphad.csv")
col1 = "Model"
col2 = "Calphad"
df["error"] = abs(df[col1] - df[col2])
df = df.sort_values(by = "error", ascending = False)

query = df[(df['error'].isnull())]
TT, FT, TF, FF = [], [], [], []
for i in df.iterrows():
	if np.isnan(i[1][col1]) and np.isnan(i[1][col2]):
		FF.append(i[1]["Alloy"])
	elif np.isnan(i[1][col1]) and not np.isnan(i[1][col2]):
		FT.append(i[1]["Alloy"])
	elif not np.isnan(i[1][col1]) and np.isnan(i[1][col2]):
		TF.append(i[1]["Alloy"])
	else:
		TT.append(i[1]["Alloy"])

# TF, TT, FF, FT = 3,18,22,2

# df_tt = pd.DataFrame.from_dict()
print({"TT": len(TT), "FT": len(FT), "TF": len(TF), "FF": len(FF)})

conf_matrix = np.array([[len(TT), len(FT)], [len(TF), len(FF)]])
print(conf_matrix)

sns.heatmap(conf_matrix, annot=True, cmap="YlGn",cbar=False)
plt.ylabel(col2)
plt.xlabel(col1)
plt.xticks([0.5, 1.5],labels = ["Miscible", "Immiscible"])
plt.yticks([0.5, 1.5], labels = ["Miscible", "Immiscible"])
plt.subplots_adjust(right=0.6, bottom=0.2)
precision = conf_matrix[0][0]/(conf_matrix[0][0]+conf_matrix[1][0])
recall = conf_matrix[0][0]/(conf_matrix[0][0]+conf_matrix[0][1])
accuracy = (conf_matrix[0][0]+conf_matrix[1][1])/np.sum(conf_matrix)
specificity = conf_matrix[1][1]/(conf_matrix[1][1]+conf_matrix[1][0])
f1 = 2*precision*recall/(precision+recall)
print(precision, recall, accuracy, specificity, f1)
plt.text(x = 2.2, y = 1.5, s = f"Precision: {precision:.2f}\n"
		 f"Recall: {recall:.2f}\n"
		 f"Accuracy: {accuracy:.2f}\n"
		 f"Specificity: {specificity:.2f}\n"
		 f"F1_score: {f1:.2f}")
plt.savefig(f"confusion_matrix_2_{col1}_{col2}_high.png", dpi = 300)

fig, ax = plt.subplots(2, 1, figsize=(8, 12))
ax[0].scatter(df["Calphad"], df["Model"], label="Misc_T")
ax[0].plot(np.linspace(200, 3201, 50), np.linspace(200, 3201, 50),)
ax[0].set_xlabel("Calphad")
ax[0].set_ylabel("Model")
ax[0].set_xlim([200, 3200])
ax[0].set_ylim([200, 3200])

# plt.subplots_adjust(hspace=0.2)
# error = abs(df["Model"] - df["Calphad"])
# sns.histplot(error, ax=ax[1], label="Error", kde=True)
# ax[1].set_xlabel("Error in T (K)")
# plt.savefig("./error_plot_high.png", dpi = 300)
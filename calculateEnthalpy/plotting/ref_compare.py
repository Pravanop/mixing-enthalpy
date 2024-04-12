
import pandas as pd

from sklearn.metrics import r2_score, mean_absolute_error


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

font = {
		'family' : 'Helvetica' ,
		'size'   : 16
		}
matplotlib.rc('font' , **font)

df = pd.read_csv("../data/input_data/ternary_tan.csv")
start = -180
stop = 180
distance = 30
x = df["Our Work"]*1000
y = df["DFT (Tan et al.)"]*0.01036*1000
x_lim = np.linspace(start, stop, 20)
x_lower = np.linspace(start + distance, stop, 20)
x_higher = np.linspace(stop - distance, start, 20)
y_higher = np.linspace(stop, start + distance, 20)
y_lower = np.linspace(start, stop - distance, 20)
plt.scatter(x,y, marker = "o", s = 80, color = "#EE7733", edgecolor = "black", linewidth = 1.2)
plt.plot(x_lim, x_lim, color = 'black', linewidth = 2, linestyle = "-")
plt.plot(x_lower, y_lower, color = 'black', linewidth = 1, linestyle = "--",
         zorder = 0, alpha = 0.5)
plt.plot(x_higher, y_higher, color = 'black', linewidth = 1, linestyle = "--",
         zorder = 0, alpha = 0.5)
# plt.fill_between(x, y - 20, y + 20, alpha=0.2)
plt.axvline(x = 0 , linewidth = 1.2 , linestyle = '-' , color = 'black' , zorder = 0, alpha = 0.8)
plt.axhline(y = 0 , linewidth = 1.2 , linestyle = '-' , color = 'black' , zorder = 0, alpha = 0.8)
errors = df.loc[abs(x-y) > distance + 15]
print(errors)
for idx, rows  in errors.iterrows():
	print(rows)
	plt.annotate(rows["Comp"], (rows["Our Work"]*1000, rows["DFT (Tan et al.)"]*0.010136*1000))
r2 = np.round(r2_score(x, y),2)
mae = np.round(mean_absolute_error(x,y),2)
plt.text(x = -75, y = 60, s = "Mean Absolute Error: \n" + str(mae) + " meV/atom" )
plt.xlim([start, stop])
plt.ylim([start, stop])
plt.xlabel("Our Work (meV/atom)")
plt.ylabel("DFT [Tan et. al.] (meV/atom)")
plt.show()
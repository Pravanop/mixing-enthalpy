import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# real = [
# 1120, 1650, 0, 1950, 475, 700, 1500, 1350, 1063, 0, 1183, 1000, 0, 1900, 1000, 1200, 900, 0,
# ]
real = [1120.0, 1650.0, 0.0, 1950.0, 475.0, 700.0, 1500.0, 1350.0, 1063.0, 0.0, 1183.0, 1000.0, 1700.0, 0.0, 1900.0, 1000.0, 1200.0, 900.0, 0.0]
# pred = [
# 1000, 1800, 0, 2200, 600, 600, 1000, 1200, 1000, 0, 800, 800, 200, 1600, 800, 200, 800, 0,
# ]
pred = [1000.0, 1700.0, 0.0, 2200.0, 600.0, 500.0, 900.0, 1200.0, 900.0, 0.0, 700.0, 800.0, 2700.0, 100.0, 1600.0, 700.0, 100.0, 800.0, 300.0]
err = abs(np.array(real)-np.array(pred))
print(len(real), len(pred))
plt.plot([min(min(real), min(pred)), max(max(real), max(pred))], [min(min(real), min(pred)), max(max(real), max(pred))], linestyle='--', zorder = 0, c= 'black')
plt.plot([min(min(real), min(pred)) + 200, max(max(real), max(pred))], [min(min(real), min(pred)), max(max(real), max(pred))-200], linestyle='--', zorder = 0, c = 'black')
plt.plot([min(min(real), min(pred)) - 200, max(max(real), max(pred))], [min(min(real), min(pred)), max(max(real), max(pred))+200], linestyle='--', zorder = 0, c= 'black')
ps = plt.scatter(real, pred, c = err, cmap = 'coolwarm', zorder = 1 )
# plt.text(s = 'Cr-W', x = 1700, y = 2200)
# plt.text(s = 'Ti-W', x = 1250, y = 200)
# plt.text(s = 'Mn-Ni', x = 1230, y = 800)
# plt.text(s = 'Ta-Zr', x =1960, y = 1600)
# plt.text(s = 'Fe-V', x =1550, y = 1000)
plt.colorbar(ps, label = 'Absolute Error (K)')
plt.ylabel('Model Prediction (K)')
plt.xlabel('Experimental Value (K)')
plt.title('Miscibility Temperature Predictions')
plt.xlim(min(min(real), min(pred)), max(max(real), max(pred)) + 100)
plt.ylim(min(min(real), min(pred)), max(max(real), max(pred)) + 100)
plt.savefig('misc_T_predictions_binary_equi.png')
plt.show()
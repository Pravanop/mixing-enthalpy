import matplotlib.pyplot as plt

import matplotlib
from IPython.core.pylabtools import figsize

matplotlib.rcParams['font.size'] = 20

x = [0.1, 0.2, 0.7]
y = [0.05, 0.15, 0.65, 0.15, 0.05, 0.05]
# 3 shades of green for colors
# 4 shades of light blue for colors
fig, ax = plt.subplots(1, 2, figsize = (12, 3))
blue_colors  = ['#00FFFF', '#00CCCC', '#009999']
red_colors = ['#FF0000', '#CC0000', '#990000']


ax[0].pie(x, labels=['$\Omega_{AB}$.$x_A$.$x_B$', '$\Omega_{BC}$.$x_B$.$x_C$', '$\Omega_{AC}$.$x_A$.$x_C$'], colors=['#00FF00', '#00CC00', '#009900'])
ax[1].pie(y, labels=['$\Omega_{AB}$.$y_A$.$y_B$', '$\Omega_{BC}$.$y_B$.$y_C$', '$\Omega_{AC}$.$y_A$.$y_C$',
					 '$\Omega_{AM}$.$y_A$.$y_M$', '$\Omega_{BM}$.$y_B$.$y_M$', '$\Omega_{CM}$.$y_C$.$y_M$'], colors=blue_colors+red_colors)
plt.subplots_adjust(wspace = 0.5)
plt.savefig('pie_chart.png')


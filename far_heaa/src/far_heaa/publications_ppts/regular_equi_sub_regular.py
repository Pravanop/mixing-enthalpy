import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# n = 2
c = ["#009988", "#EE7733", "#0077BB", "#CC3311"]
legend = ['Binary', 'Ternary', 'Quartenary', 'Quinary']
fig, ax = plt.subplots(1, 4, sharex=True, sharey=True)
for idx, n in enumerate([2, 3, 4, 5]):
    regular = pd.read_csv(f'../output_data/predictions/regular/misc_T_{n}_BCC.csv')
    sub_regular = pd.read_csv(f'../output_data/predictions/sub-regular/misc_T_{n}_BCC.csv')


    y_Reg = regular['Misc_Temp'].to_numpy()
    y_subReg = sub_regular['Misc_Temp'].to_numpy()
    sns.histplot(abs(y_Reg-y_subReg), ax = ax[idx], color = c[idx], kde = True)
    ax[idx].set_xlim([0, 750])
    ax[idx].set_title(legend[idx])
    # plt.scatter(y_Reg, y_subReg, c = c[idx], zorder =1)


fig.supxlabel('MAE (K)')
plt.subplots_adjust(wspace = 0, hspace=0)
plt.show()
# plt.xlabel('Regular Model')
# plt.ylabel('Sub Regular Model')
# plt.title('Miscible Temperature Using Different Models')
# plt.plot(range(0, 3000), range(0, 3000), c = 'black',linestyle = '--', zorder = 0)
# plt.xlim([0, 3000])
# plt.ylim([0, 3000])
# plt.legend(['Binary', 'Ternary', 'Quartenary', 'Quinary'])
plt.show()
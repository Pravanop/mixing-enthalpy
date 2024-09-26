import numpy as np
import pandas as pd
import seaborn as sns
from far_heaa.high_throughput.equi_miscT_predictions import EquiMiscTPredictions
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.math_operations.thermo_calculations import ThermoMaths
import matplotlib.pyplot as plt
df = pd.read_csv('/Users/mcube/Desktop/Projects/calculateEnthalpy/calphad/expt_binary_prediction.csv')
print(df.head())
compositions = df['0']
expt_preds = df['expt_T_misc']
print(len(expt_preds))
mH = MetadataHandler()
meta_data = mH.get_metadata
meta_data['flags']['equi_flag'] = False
meta_data['flags']['corrections'] = False
# meta_data['flags']['im_flag'] = False
equi_list = []
for i in compositions:
    try:
        composition = i.split('-')
    except:
        continue
    equi_miscT = EquiMiscTPredictions(dim = 2,
                                      element_list = composition,
                                      lattice='min',
                                      meta_data=meta_data,
                                      save_flag = False,
                                      generate_plot=False)

    answer = float(equi_miscT.make_predictions()['Misc_Temp'])
    avg_tm = ThermoMaths().avg_T_melt(composition=composition, mol_ratio= [0.5, 0.5])
    if answer == avg_tm:
        equi_list.append(np.nan)
    else:
        equi_list.append(answer)

confusion_matrix = np.zeros((2,2))

pred_filt = []
expt_filt = []
for idx, i in enumerate(equi_list):
    true = np.isnan(expt_preds[idx])
    pred = np.isnan(i)
    if true and pred:
        confusion_matrix[1][1] += 1
    if true and not pred:
        confusion_matrix[1][0] += 1
    if not true and pred:
        confusion_matrix[0][1] += 1
    if not true and not pred:
        pred_filt.append(i)
        expt_filt.append(expt_preds[idx])
        confusion_matrix[0][0] += 1

print(pred_filt, expt_filt)
# confusion_matrix = confusion_matrix.astype(int)
# sns.heatmap(confusion_matrix, annot=True, fmt='d')
# plt.show()

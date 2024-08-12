import numpy as np
import pandas as pd
import xlsxwriter
wb = xlsxwriter.Workbook('4element.xlsx',
                         {'strings_to_numbers': True})
sheets = ["BCC", "FCC"]
for n in range(2):
    sheet = pd.read_excel('RHEAs_phase_prediction-main/enthalpy_data_and_predictions/bokas.xlsx',
                          sheet_name=sheets[n],
                          index_col=0,
                          usecols='A:AA',
                          nrows=26)
    s1 = wb.add_worksheet('Quaternary ' + sheets[n])
    s1.write(0, 0, "comp")
    for i in range(4):
        s1.write(0, i + 1, "e" + str(i+1))
    s1.write(0, 5, "enthalpy")
    s1.write(0, 6, "entropy")
    cnt = 1
    for e1 in range(26):
        a = sheet.columns.values[e1]
        for e2 in range(e1 + 1, 26):
            b = sheet.columns.values[e2]
            for e3 in range(e2 + 1, 26):
                c = sheet.columns.values[e3]
                for e4 in range(e3 + 1, 26):
                    d = sheet.columns.values[e4]
                    s1.write(cnt, 1, a)
                    s1.write(cnt, 2, b)
                    s1.write(cnt, 3, c)
                    s1.write(cnt, 4, d)
                    s1.write(cnt, 0, str(a + b + c + d))
                    sm = 0
                    var = [e1, e2, e3, e4]
                    for i in range(4):
                        for j in range(i + 1, 4):
                            add = sheet.iloc[var[j], var[i]]
                            sm += float(add)
                    s1.write(cnt, 5, sm/16*4)
                    s1.write(cnt, 6, 0.000119)
                    cnt += 1

wb.close()

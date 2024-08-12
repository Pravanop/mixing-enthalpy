import numpy as np
import pandas as pd
import xlsxwriter
wb = xlsxwriter.Workbook('3element.xlsx',
                         {'strings_to_numbers': True})
sheets = ["BCC", "FCC"]
for n in range(2):
    sheet = pd.read_excel('RHEAs_phase_prediction-main/bokas.xlsx',
                          sheet_name=sheets[n],
                          index_col=0,
                          usecols='A:AA',
                          nrows=26)
    s1 = wb.add_worksheet('Ternary ' + sheets[n])
    s1.write(0, 0, "comp")
    for i in range(4):
        s1.write(0, i + 1, "e" + str(i+1))
    s1.write(0, 4, "enthalpy")
    s1.write(0, 5, "entropy")
    cnt = 1
    for e1 in range(26):
        a = sheet.columns.values[e1]
        for e2 in range(e1 + 1, 26):
            b = sheet.columns.values[e2]
            for e3 in range(e2 + 1, 26):
                c = sheet.columns.values[e3]
                s1.write(cnt, 1, a)
                s1.write(cnt, 2, b)
                s1.write(cnt, 3, c)
                s1.write(cnt, 0, str(a + b + c))
                sm = 0
                var = [e1, e2, e3]
                for i in range(len(var)):
                    for j in range(i + 1, len(var)):
                        add = sheet.iloc[var[j], var[i]]
                        sm += float(add)
                s1.write(cnt, 4, sm/9*4) # simplify into common math enthalpy function
                s1.write(cnt, 5, 0.0000947)
                cnt += 1

wb.close()

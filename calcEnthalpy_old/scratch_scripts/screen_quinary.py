import numpy as np
import pandas as pd
import xlsxwriter
sheets = ["BCC", "FCC"]
wb = xlsxwriter.Workbook('5z.xlsx',
                         {'strings_to_numbers': True})
for n in range(1):
    sheet = pd.read_excel('zhaohan/datasheet.xlsx',
                          index_col=0)
    # s1 = wb.add_worksheet("Quinary " + sheets[n])
    length = len(sheet.index)
    s1 = wb.add_worksheet()
    s1.write(0, 0, "comp")
    for i in range(5):
        s1.write(0, i + 1, "e" + str(i+1))
    s1.write(0, 6, "enthalpy")
    s1.write(0, 7, "entropy")
    cnt = 1
    for e1 in range(length):
        a = sheet.columns.values[e1]
        for e2 in range(e1 + 1, length):
            b = sheet.columns.values[e2]
            for e3 in range(e2 + 1, length):
                c = sheet.columns.values[e3]
                for e4 in range(e3 + 1, length):
                    d = sheet.columns.values[e4]
                    for e5 in range(e4 + 1, length):
                        e = sheet.columns.values[e5]
                        s1.write(cnt, 1, a)
                        s1.write(cnt, 2, b)
                        s1.write(cnt, 3, c)
                        s1.write(cnt, 4, d)
                        s1.write(cnt, 5, e)
                        s1.write(cnt, 0, str(a + b + c + d + e))
                        sm = 0
                        var = [e1, e2, e3, e4, e5]
                        for i in range(5):
                            for j in range(i + 1, 5):
                                add = sheet.iloc[var[j], var[i]]
                                sm += float(add)
                        s1.write(cnt, 6, sm/25*4)
                        s1.write(cnt, 7, 0.000139)
                        cnt += 1

wb.close()

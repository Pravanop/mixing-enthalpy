def list_to_pairs(elements:list):
    out = []
    for i in elements:
        for j in elements:
            if i is not j:
                out.append((i, j))
    return out

if __name__ == '__main__':
    els = ['Al', 'Co', 'Cr', 'Cu', 'Fe', 'Hf', 'Mn', 'Mo',
           'Nb', 'Ni', 'Ta', 'Ti', 'W', 'Zr', 'V', 'Mg', 'Re',
           'Os', 'Rh', 'Ir', 'Pd', 'Pt', 'Ag', 'Au', 'Zn', 'Cd', 'Ru']
    print(list_to_pairs(els))
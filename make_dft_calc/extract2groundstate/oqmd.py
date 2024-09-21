import qmpy_rester as qr

## Return list of data
with qr.QMPYRester() as q:
    kwargs = {
        'element_set': '(Fe-Mn),O',      # composition include (Fe OR Mn) AND O
        'stability': '<-0.1',            # hull distance smaller than -0.1 eV
        'natom': '<10',                  # number of atoms less than 10
        }
    list_of_data = q.get_oqmd_phases(**kwargs)
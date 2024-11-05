import numpy as np
from matplotlib import pyplot as plt
import mpltern
from sympy import symbols
from sympy.functions.elementary.exponential import log
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices.dense import Matrix
from far_heaa.grids_and_combinations.combination_generation import MultinaryCombinations
from far_heaa.grids_and_combinations.grid_creation import CompositionGrid
from far_heaa.math_operations.thermo_calculations import ThermoMaths
from far_heaa.io.metadata_handler import MetadataHandler
from far_heaa.visualization_toolkit.visualizations import Visualizations

tm = ThermoMaths()
def create_xsyms(n):
    '''returns a string of symbols: (x1,x2,...,xn-1)
    Contributor: John Cavin'''
    s = ''
    for i in range(n-1):
        s=s+'x{} '.format(int(i+1))
    return symbols(s)

def x_N(syms):
    '''input: n-1 syms
    returns 1-x1-x2-...x_(n-1)
    Contributor: John Cavin'''
    xn = parse_expr('1')
    #print(syms)
    for s in syms:
        xn = xn-s
    return xn


def create_Ssym(xsyms):
    out = parse_expr('0')
    for i in xsyms:
        out += - tm.kb * i * log(i + 1e-4)

    return out

def create_Hsym(xsyms, composition, data, lattice):
    binaries = list(MultinaryCombinations.create_multinary(composition, no_comb = [2]).values())[0]
    out = parse_expr('0')
    for binary in binaries:
       omega = data[binary][lattice]
       two_eles = binary.split('-')
       out += omega*xsyms[composition.index(two_eles[0])]*xsyms[composition.index(two_eles[1])]

    return out

def hessian(f,syms):
    '''compute the hessian of a matrix of symbols
    Contributor: John Cavin'''
    out = Matrix([[f.diff(x).diff(y) for x in syms] for y in syms])
    return out


def G_sym(H_sym, S_sym, T):
    return H_sym - T*S_sym

def find_eigenvalue(H_num):
    H_num = np.array(H_num).astype(np.float64)
    w, _ = np.linalg.eig(H_num)
    return w

mH = MetadataHandler()
meta_data = mH.get_metadata
viz = Visualizations(lattice= "BCC", meta_data = meta_data)
composition = ['Cr', 'V', 'Ti']


xsym = list(create_xsyms(len(composition)))
xsyms = xsym + [x_N(xsym)]
H_sym = create_Hsym(xsyms, composition, viz.data, 'BCC')
S_sym = create_Ssym(xsyms)


T = symbols("T")


G = G_sym(H_sym, S_sym, T)
hess = hessian(G, xsym)

# mol = [0, 0.5]
# temp_eigen = []
# temp_grid = np.linspace(100, 2800, 20)
# for temperature in temp_grid:
#     H_num = hess.subs(list(zip(xsym + [T], mol + [temperature])))
#     eigen = find_eigenvalue(H_num)
#     temp_eigen.append(eigen)
#     if np.all(eigen > 0):
#         print(temperature)
#
# temp_eigen = np.array(temp_eigen)
# # for i in range(temp_eigen.shape[1]):
# #     plt.plot(temp_grid, temp_eigen[:, i], marker = 'o')
#
#
# plt.plot(temp_grid, temp_eigen[:, 0], marker = 'o')
# plt.axhline(y = 0)
# plt.show()

temperature = 100
mol_grid = CompositionGrid.create_mol_grid(len(composition), 40)

unstable = []
meta_stable = []
for mol in mol_grid:
    mol = list(mol[:len(composition)-1])
    H_num = hess.subs(list(zip(xsym + [T], mol + [temperature])))
    eigen = find_eigenvalue(H_num)
    is_metastable = True
    is_unstable = False
    for i in eigen:
        if i > 0:
            continue
        if i < 0:
            is_unstable = True
            is_metastable = False

    if is_metastable:
        meta_stable.append(mol)
    if is_unstable:
        unstable.append(mol)

fig = plt.figure()

ax = fig.add_subplot(projection="ternary")
ax.grid()

for i in unstable:
    ax.scatter(i[0], i[1],1 - i[0] - i[1], c = 'red')

for i in meta_stable:
    ax.scatter(i[0],  i[1],  1 - i[0] - i[1], c = 'blue')

plt.show()








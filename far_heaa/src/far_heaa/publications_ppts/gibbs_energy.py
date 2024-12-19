

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size': 12})
x = np.linspace(0.0, 1, 100)


h = 0.1*x*(1-x)



plt.plot(x, h, c = 'black')
#shades of blue
gibbs_colors = ['blue', 'skyblue', 'deepskyblue']
#shades of red
s_colors = ['red', 'salmon', 'darkred']
idx = 0
for T in [200, 400, 500]:
	
	s = T*8.315e-5*x*np.log(x+1e-9) + T*8.315e-5*(1-x)*np.log(1-x + 1e-9)
	# plt.plot(x, s, c = s_colors[idx])
	plt.plot(x, h + s, c = gibbs_colors[idx], linewidth = 2)
	idx += 1

plt.title('Gibbs Energy with Increasing Temperature')
plt.xlabel('x')
plt.ylabel('Gibbs Energy (ev/atom)')
plt.axhline(0, c='black', label = "_no_legend_")

plt.legend([0] + [200, 400, 500])
plt.savefig('gibbs_energy.png')
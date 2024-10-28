import seaborn as sns
import matplotlib.pyplot as plt

x = ['Total Paths', 'Miscible to Miscible', 'Immiscible to Immiscible', 'Miscible to Immiscible', 'Immiscible to Miscible']
y = [840, 417, 217, 108, 98]

sns.barplot(x= y, y = x, color = '#BB5566', alpha = 1)
plt.yticks(x)
plt.xlabel('# Paths')
plt.subplots_adjust(left = 0.3, right = 0.95)
plt.show()
from partitioned import AlloySystem as alloy
import pandas as pd

newalloy = alloy(directory="FCC")
newalloy.run_ehull_combinations(5, 1000)
# libraries
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd
import math
import operator
import itertools
import os.path
import sys
import ast

import APQ_Simulation_Engine as APQ

# full simulation, takes about 4it's .5 hours to run

num_customers = 6000
num_sims = 50
lamvals = [(0.5,0.3), (0.2,0.7)]
bdvals = [(b,d) for b in [0.2,0.4,0.6,0.8] for d in [0,2,6]] + [(0,0)] + [(1,0)]

np.random.seed(4811)
APQ.multiapq(num_customers, num_sims, lamvals, bdvals)


# simulate a two level M/M/1 delayed accumulating priority queue

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
import time
from tqdm import tqdm, tqdm_gui
from timeit import default_timer as timer

# helper lemma
def cummean(x):  
    return(np.cumsum(x)/np.arange(1, len(x)+1))

# from customer dataframe get customer with most priority
# only works if priority preserved within classes
def priority(df, time, b, d):
    df1 = df[df.level==1]
    df2 = df[df.level==2]
    if (len(df1)==0 or len(df2)==0):
        return(min(df.index))
    c1 = min(df1.index)
    c2 = min(df2.index)
    c1p = (time - df.loc[c1,'arrival'] - 0) * 1
    c2p = (time - df.loc[c2,'arrival'] - d) * b
    if c1p >= c2p:
        return(c1)
    return(c2)

##################
# run simulation with multiple types of queueing disciplines and runs
##################

def multiapq(num_customers, num_sims, lamvals, bdvals):

    # path to export files to
    basepath = '/Users/blairbilodeau/Documents/Research/Projects/Delayed APQ/delayed-apq/code/APQ_Simulation_Results/'
    sys.path.append(basepath)

    # for easier param reference
    mu = 1 # for simplicity

    for sim in range(num_sims):

        for (lam1,lam2) in lamvals: 

            lam = lam1 + lam2
            class2prob = lam2/lam
        
            # compare apples to apples by using same arrival data for each queueing type
            level = np.random.binomial(1, class2prob, num_customers) + 1
            interarrival = np.random.exponential(1/lam, num_customers)
            service_length = np.random.exponential(1/mu, num_customers)
            arrival = np.array(np.cumsum(interarrival)) 

            for (b,d) in bdvals:

                # export file
                filename = 'APQsim' + '_ncust_' + str(num_customers) + '__' + 'lam1_' + str(lam1) + '_lam2_' + str(lam2) + '__' + 'b_' + str(b) + '_d_' + str(d) + '__' + 'r_' + str(sim) + '.csv'
                filepath = basepath + filename

                if os.path.isfile(filepath):
                    print(filename + ' already exists')
                else:

                    # start with empty customer dataframe, idle server, and time 0
                    customers = pd.DataFrame(index=range(num_customers), columns=['level', 'interarrival', 'service_length', 'arrival', 'service', 'departure', 'wait'])

                    # fill in inter-arrivals, service lengths, and class types
                    customers.level = level
                    customers.interarrival = interarrival
                    customers.service_length = service_length
                    customers.arrival = arrival

                    # move to first departure
                    customers.loc[0, 'service'] = customers.loc[0, 'arrival']
                    time = customers.loc[0, 'service'] + customers.loc[0, 'service_length']
                    customers.loc[0, 'departure'] = time 
                    min_cust = 0
                    last_dep = 0

                    # process each departure except the last one
                    for departures in range(num_customers-2):

                        # find the customer index of the most recent departure and the arrival time of the customer after them
                        max_cust = min(last_dep, num_customers-2)
                        max_arrival = customers.loc[max_cust+1, 'arrival']
                        while max_arrival <= time and max_cust < (num_customers-2):
                            max_cust += 1
                            max_arrival = customers.loc[max_cust+1, 'arrival']

                        # proxy for customers who have arrived
                        arrived_customers = customers[min_cust:(max_cust+1)]

                        # remove those who have not yet entered service
                        if len(arrived_customers) > 0:
                            active_customers = arrived_customers[arrived_customers.service != arrived_customers.service]
                        else:
                            active_customers = arrived_customers

                        # if the queue is empty
                        if len(active_customers)==0:
                            # next customer to arrive is right after the most recent arrival
                            cust = max_cust + 1
                            customers.loc[cust, 'service'] = customers.loc[cust, 'arrival']
                            min_cust = cust # the next person can't have arrived before them since they found the queue empty

                        # if somebody eligible to enter service
                        else:
                            # customer with highest priority
                            cust = priority(active_customers, time, b, d)
                            customers.loc[cust, 'service'] = time
                            if customers.loc[cust,'level']==2:
                                min_cust = cust # a level 2 can't be served before someone who arrived before them
                        
                        time = customers.loc[cust, 'service'] + customers.loc[cust, 'service_length'] 
                        customers.loc[cust, 'departure'] = time
                        last_dep = cust

                    # compute averages
                    customers.wait = customers.service - customers.arrival

                    print('NCust: ' + str(num_customers) + ', lam1: ' + str(lam1) + ', lam2: ' + str(lam2) + ', b: ' + str(b) + ', d: ' + str(d) + ', run: ' + str(sim+1) + '/' + str(num_sims))
                    customers.to_csv(filepath)

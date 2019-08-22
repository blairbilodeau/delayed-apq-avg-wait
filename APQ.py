# Class for two level delayed accumulating priority queue
# Blair Bilodeau
# simulate an accumulating priority queue (with delays)

# libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import math
import operator
import os.path
import sys
import ast
import time
from tqdm import tqdm, tqdm_gui
from timeit import default_timer as timer

# parameters
mu = 1
lam1 = 0.5
lam2 = 0.3
b = [1, 0] # rate at which customers accumulate
d = [0, 0] # initial delay before customers accumulate

lam = lam1 + lam2
rho1 = lam1 / mu
rho2 = lam2 / mu
rho = rho1 + rho2
# lam1A = lam1 * (1-b[1]/b[0])

# from customer dataframe get customer with most priority
def priority(df, time):
    accum = (time - df.arrival - np.array([d[level - 1] for level in df.level])) * np.array([b[level - 1] for level in df.level])
    ind = np.argmax(accum)
    return(ind)

cummean = lambda x:  np.cumsum(x)/np.arange(1, len(x)+1)


# start with empty customer dataframe, idle server, and time 0
customers = pd.DataFrame(index=range(num_customers), columns=['level', 'interarrival', 'service_length', 'arrival', 'service', 'departure', 'wait'])

# fill in inter-arrivals, service lengths, and class types
customers.level = np.random.binomial(1,rho2,num_customers) + 1
customers.interarrival = np.random.exponential(1/lam, num_customers)
customers.service_length = np.random.exponential(1/mu, num_customers)
arrival = np.array(np.cumsum(customers.interarrival))
customers.arrival = arrival

# move to first departure
customers.loc[0, 'service'] = customers.loc[0, 'arrival']
time = customers.loc[0, 'service'] + customers.loc[0, 'service_length']
customers.loc[0, 'departure'] = time 


time_df = pd.DataFrame(index=range(num_customers-1), columns=['t1', 't2', 't3', 't4'])

start = timer()
# process each departure except the last one
for departure in range(num_customers-1):
    # customers who have already arrived and not yet entered service
    start1 = timer() # this is the slow part!!
    #active = customers.query('arrival < @time and service != service') # nan is never equal to itself
    active = customers[np.all([list(customers.arrival < time),list(customers.service != customers.service)], axis=0)]
    end1 = timer()
    time_df.loc[departure, 't1'] = end1 - start1

    # if nobody eligible to enter service
    if len(active)==0:
        # next customer to arrive
        start4 = timer()
        cust = min(customers.query('service != service').index)
        customers.set_value(cust, 'service', customers.values[cust, 3]) # arrival
        time = customers.values[cust, 4] + customers.values[cust, 2] # service + service_length
        customers.set_value(cust, 'departure', time)
        end4 = timer()
        time_df.loc[departure, 't4'] = end4 - start4
        time_df.loc[departure, 't2'] = 0
        time_df.loc[departure, 't3'] = 0

    # if somebody eligible to enter service
    else:
        # customer with highest priority
        start2 = timer()
        cust = priority(active, time)
        end2 = timer()
        start3 = timer()
        customers.set_value(cust, 'service', time)
        time = customers.values[cust, 4] + customers.values[cust, 2] # service + service_length
        customers.set_value(cust, 'departure', time)
        # customer = customers.loc[cust].copy()
        # customer.service = time
        # time = customer.service + customer.service_length
        # customer.departure = time 
        # customers.loc[cust] = customer
        end3 = timer()
        time_df.loc[departure, 't3'] = end3 - start3
        time_df.loc[departure, 't2'] = end2 - start2
        time_df.loc[departure, 't4'] = 0

customers.wait = customers.service - customers.arrival
end = timer()
total = end - start

mean1 = np.mean(time_df[time_df.t1 != 0].t1)
mean2 = np.mean(time_df[time_df.t2 != 0].t2)
mean3 = np.mean(time_df[time_df.t3 != 0].t3)
mean4 = np.mean(time_df[time_df.t4 != 0].t4)
total_sum = sum(time_df.t1) + sum(time_df.t2) + sum(time_df.t3) + sum(time_df.t4)



plt.plot(customers.index, customers.wait, customers.index, cummean(customers.wait), 'r--')
plt.show()

def apq(num_customers, burn_in, burn_out, num_reps):

    avg_wait = pd.DataFrame(index=range(num_reps), columns=['level1', 'level2', 'overall'])

    for rep in range(int(num_reps)):
        # start with empty customer dataframe, idle server, and time 0
        customers = pd.DataFrame(index=range(num_customers), columns=['level', 'interarrival', 'service_length', 'arrival', 'service', 'departure', 'wait'])

        # fill in inter-arrivals, service lengths, and class types
        customers.level = np.random.binomial(1,rho2,num_customers) + 1
        customers.interarrival = np.random.exponential(1/lam, num_customers)
        #customers.service_length = np.random.exponential(1/mu, num_customers)
        customers.service_length = np.repeat(1/mu, num_customers)
        arrival = np.array(np.cumsum(customers.interarrival))
        customers.arrival = arrival

        # move to first departure
        customers.loc[0, 'service'] = customers.loc[0, 'arrival']
        time = customers.loc[0, 'service'] + customers.loc[0, 'service_length']
        customers.loc[0, 'departure'] = time 

        # process each departure except the last one
        for departure in range(int(num_customers-1)):

            if departure % 100 == 0:
                print('rep:' + str(rep) + ', dep:' + str(departure))

            # customers who have already arrived and not yet entered service
            active = customers[np.all([list(customers.arrival < time),list(customers.service != customers.service)], axis=0)]

            # if nobody eligible to enter service
            if len(active)==0:
                # next customer to arrive
                cust = min(customers.query('service != service').index)
                customers.set_value(cust, 'service', customers.values[cust, 3]) # arrival
                time = customers.values[cust, 4] + customers.values[cust, 2] # service + service_length
                customers.set_value(cust, 'departure', time)

            # if somebody eligible to enter service
            else:
                # customer with highest priority
                cust = priority(active, time)
                customers.set_value(cust, 'service', time)
                time = customers.values[cust, 4] + customers.values[cust, 2] # service + service_length
                customers.set_value(cust, 'departure', time)

        # compute averages
        customers.wait = customers.service - customers.arrival
        customers = customers[burn_in:(num_customers-burn_out)]
        avg_wait.loc[rep, ['level1', 'level2', 'overall']] = [np.mean(customers[customers.level==1].wait), np.mean(customers[customers.level==2].wait), np.mean(customers.wait)]

    return([customers, avg_wait])

test = apq(10000,750,250,1)

test_customers = test[0]
test_wait = test[1]
m1 = np.mean(test_wait.level1)
m2 = np.mean(test_wait.level2)
m = np.mean(test_wait.overall)

plt.plot(test_customers.index, test_customers.wait, test_customers.index, cummean(test_customers.wait), 'r--')
plt.show()

plt.plot(test_customers[test_customers.level==1].index, test_customers[test_customers.level==1].wait, test_customers[test_customers.level==1].index, cummean(test_customers[test_customers.level==1].wait), 'r--')
plt.show()

plt.plot(test_customers[test_customers.level==2].index, test_customers[test_customers.level==2].wait, test_customers[test_customers.level==2].index, cummean(test_customers[test_customers.level==2].wait), 'r--')
plt.show()


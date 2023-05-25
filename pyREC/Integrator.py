import numpy as np
import pandas as pd


def integrator(dataseries, unit, period):
    """

    :param dataseries: full simulation for one year
    :param unit: 'energy' or 'power'
    :param period: analysis period ('year','month','day')
    :return: integral over analysis period
    """
    time=8760/len(dataseries)
    dataseries = np.array(dataseries)
    r_t = int(1 / time)
    month = [0, 744 * r_t, 1416 * 1 * r_t, 2160 * r_t, 2880 * r_t, 3624 * r_t, 4344 * r_t, 5088 * r_t, 5832 * r_t,
             6552 * r_t,
             7296 * r_t, 8016 * r_t, 8760 * r_t]

    if period == 'year':
        if unit == 'energy':
            tot_period = dataseries.sum()
        else:
            for i in range(len(dataseries)):
                dataseries[i] = dataseries[i] * time
            tot_period = dataseries.sum()

    elif period == 'month':
        tot_period = []
        for j in range(len(month) - 1):
            tot_period_i = 0
            for i in range(month[j], month[j + 1]):
                if unit == 'power':
                    tot_period_i += dataseries[i] * time
                else:
                    tot_period_i += dataseries[i]
            tot_period.append(tot_period_i)

    elif period == 'day':
        tot_period = []
        for j in range(int(len(dataseries) / (24 * r_t))):
            tot_period_i = 0
            for i in range(j * 24 * r_t, j * 24 * r_t + 24 * r_t):
                if unit == 'power':
                    tot_period_i += dataseries[i] * time
                else:
                    tot_period_i += dataseries[i]

            tot_period.append(tot_period_i)
    else:
        tot_period=None

    tot_period=np.array(tot_period)

    return tot_period



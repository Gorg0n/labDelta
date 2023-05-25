import numpy as np


class User:
    def __init__(self, id, carrier,dem, pod=None, group=None, ma=None, ts=None, plants=None, space=None):
        """
        to simulate an electricity consumer
        :param id: identification code
        :param carrier: energy vector 'electricity' or 'heat'
        :param dem: load curve [kW]
        :param pod: point of delivery
        :param group: residential,commercial,industrial,other type.
        :param ma: electricity market area
        :param ts: transformer station
        :param plants: list of plants to which the user is physically connected
        :param space: spaces available for the installation of new systems [m2]
        """
        self.id = id
        self.carrier=carrier
        self.dem=np.array(dem)
        if self.carrier=='electricity':
            self.pod = pod
            self.ma = ma
            self.ts = ts

        self.group = group
        self.plants = plants
        self.space = space
        self.en_perf_evolution = {}  # {'dem':initial demand ,'deficit': residual demand}
        if len(carrier)>1:
            for i in range(len(carrier)):
                self.en_perf_evolution[carrier[i]] = np.array(dem[i])
        else:
            self.en_perf_evolution[carrier[0]] = np.array(dem)








import numpy as np
from scipy.signal import argrelextrema
from .river import River

class WillettSPM(object):
    
    def __init__(self, surface, precipitation_rate, erodibility, uplift):
        self.surface = surface
        self.precipitation_rate = precipitation_rate
        self.erodibility = erodibility
        self.uplift = uplift

        self.rivers = self._get_rivers()
        
    def _get_maximas(self):
        y = self.surface[:, 1]
        return argrelextrema(y, np.greater_equal, mode='clip', order=2)[0]
        
    def _get_minimas(self):
        y = self.surface[:, 1]
        return argrelextrema(y, np.less_equal, mode='clip', order=2)[0]
    
    def _get_rivers(self):
        rivers = []
        self.minimas = minimas = self._get_minimas()
        self.maximas = maximas = self._get_maximas()
        P = self.precipitation_rate
        K = self.erodibility
        U = self.uplift
        
        for node in maximas:
            idx = np.searchsorted(minimas, node)

            if idx - 1 >=0:
                indices = np.arange(node, minimas[idx-1] - 1, -1)
                if indices.size > 0:
                    river = River(self.surface[indices],indices,
                                  node, minimas[idx-1], P[indices], K[indices], U[indices])
                rivers.append(river)
            if idx < len(minimas):
                indices = np.arange(node, minimas[idx] + 1, 1)
                if indices.size > 0:
                    river = River(self.surface[indices],indices,
                                  node, minimas[idx], P[indices], K[indices], U[indices])
                rivers.append(river)
            if idx == len(minimas):
                indices = np.arange(node, len(minimas), 1)
                if indices.size > 0:
                    river = River(self.surface[indices], indices,
                                  node, len(minimas), P[indices], K[indices], U[indices])
                rivers.append(river)
                
        return rivers    
    
    def plot_rivers(self):
        import matplotlib.pyplot as plt
        rivers = self.rivers
        for river in rivers:
            plt.plot(river.surface[:, 0], river.surface[:, 1])

    def plot_surface(self):
        import matplotlib.pyplot as plt
        plt.plot(self.surface[:,0], self.surface[:,1])
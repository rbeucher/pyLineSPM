import numpy as np
import sys


class River(object):
    
    def __init__(self, surface, nodes, max_node, min_node, precipitation_rate, erodibility, uplift, sea_level=0.):

        self.surface = surface
        # The surface profile starts at the water divide down to the bottom node of the river.
        # We use a local reference frane where the origin is set at the water divide.
        self.x = np.abs(self.surface[:,0] - self.surface[0, 0]) # Might want to generalise this.
        self.y = self.surface[:,1]
        self.nodes = nodes
        self.max_node = max_node
        self.min_node = min_node
        self.precipitation_rate = precipitation_rate
        self.erodibility = erodibility
        self.uplift = uplift
        self.sea_level = 0.

        # Cunulative length along the river
        self.dL = self.cumlength = np.cumsum(np.diff(self.x)**2 + np.diff(self.y)**2)
        # Total cumulative length
        self.length = self.dL[-1]
        
        if np.argmin(self.y) != len(self.y) - 1:
            raise ValueError(f"Minimum elevation should be at node {len(self.y) - 1} not {np.argmin(self.y)}")
        if np.argmax(self.y) != 0:
            raise ValueError(f"Maximum elevation should be at node 0")
        if not np.allclose(np.diff(self.x), np.diff(self.x)[0]):
            raise ValueError("Spacing in x is not uniform")
        if not np.all(np.diff(self.x) >= 0.):
            raise ValueError("X coordinates should increase")
        for key, val in {"nodes": nodes, "precipitation_rate": precipitation_rate, "erodibility": erodibility, "uplift": uplift}.items():
            if val.shape != surface[:,0].shape:
                raise ValueError(f"Wrong array shape {val.shape} for array '{key}'")
        
    def get_discharge(self, ka=0.3, h=2.0):
        # The Precipitation is integrated to define discharge
        # at each node. Assuming that precipitation rate is
        # a function of x, the discharge is calculated by integrating in
        # x, from the water divide to the node of interest.
        P = self.precipitation_rate
        x = self.x
        dx = np.diff(x)[0] # Here we assume constant spacing in x.
        Psum = np.cumsum(P * dx) 
        return Psum * ka * x**(h-1)
    
    def calculate_zT(self, m, n, ka, h, theta, xC, yC, dt):
        """ Calculate and apply erosion in the trunk channel"""
        
        # The Lowest node in the river is assumed to
        # be fixed, the erosion rate is equal to 0.
        # We calculate the value of the erosion rate
        # starting from the second lowest node and moving
        # in the upstream direction.
        
        # We reverse all arrays
        K = self.erodibility[::-1]
        Q = self.get_discharge(ka, h)[::-1]
        zr = self.calculate_zR(m, n, ka, h, theta, xC, yC, dt)
        zt = self.y - zr / 2.0
        zt = zt[::-1]
        x = self.x[::-1]
        dx = np.abs(np.diff(x)[0])

        ## Get first index above sea level
        first_idx = np.argmax(zt - self.sea_level > 0.)
        print(zt[first_idx])
        
        e0 = 1e-6
        erosion_rate = np.zeros_like(zt)
        
        # We integrate from the boundary node to the
        # main divide.
        for idx in range(first_idx, len(x)):
            if x[idx] > xC:
                rtol = 1e3
                edot = e0
                while rtol > 1e-7:
                    prev = edot
                    edot = -K[idx] * Q[idx]**m * (zt[idx] - zt[idx - 1] + dt * (edot - erosion_rate[idx-1]))**n / dx**n
                    rtol = np.abs((edot - prev) / prev)
                erosion_rate[idx] = edot
            zt[idx] = zt[idx] + dt * erosion_rate[idx]

        self.zT = zt[::-1]
        return self.zT

    @staticmethod
    def calculate_yD(ka, x, h):
        return 0.5 * ka * x**(h-1.0)

    @staticmethod
    def calculate_N(m, n, ka, h, K, P, xD):
        return (1.0 / (xD**(h*m) * K * P**m * ka**m))**(1/n)

    @staticmethod
    def calculate_D(m, n, h, N, yC, yD):
        hmn = h*m/n
        D = 0.
        with np.errstate(divide='ignore', invalid='ignore'):
            if hmn == 1.0:
                D = np.where(yD == 0, 0., -0.5 * N*np.log(yC/yD))
            else:
                D = np.where(yD == 0, 0., 0.5 * yD**(1-hmn) * N * (1-hmn)**-1 * (1-(yC/yD)**(1-hmn)))
        return D

    def calculate_zR(self, m, n, ka, h, theta, xC, yC, dt):
        
        theta = np.radians(theta)
        K = self.erodibility
        P = self.precipitation_rate
        U = self.uplift

        x = self.x        
        xD = np.max(x)
        zR = np.zeros_like(x)
        
        hmn = h*m/n

        yD = self.calculate_yD(ka, x, h)
        N = self.calculate_N(m, n, ka, h, K, P, xD)
        D = self.calculate_D(m, n, h, N, yC, yD)
        zR = np.where(yD > yC, yC * np.tan(theta) + 2. * D * U**(1/n)*xD**hmn, yD * np.tan(theta))

        return zR
    
    def calculate_zD(self, m, n, ka, h, theta, xC, yC, dt):
        
        theta = np.radians(theta)
        K = self.erodibility
        P = self.precipitation_rate
        U = self.uplift
        
        x = self.x
        xD = np.max(x)
        zT = self.zT
        zD = np.zeros_like(self.zT)
        
        hmn = h*m/n

        yD = self.calculate_yD(ka, x, h)
        N = self.calculate_N(m, n, ka, h, K, P, xD)
        D = self.calculate_D(m, n, h, N, yC, yD)

        with np.errstate(divide='ignore', invalid='ignore'):
            zD = np.where(yD > yC, zT + yC * np.tan(theta) + 2. * D * U**(1/n)*xD**hmn, zT + yD * np.tan(theta))
            edotR = np.where(D, -K * P**m*ka**m * ((zD - zT - yC * np.tan(theta)) / (2*D*xD**hmn))**n, 0)
        zD = zD + dt * edotR

        return zD
    
    def update_surface(self, m=0.5, n=1.0, ka=0.3, h=2.0, theta=30, xC=10, yC=10, dt=0.1):
        zT = self.calculate_zT(m, n, ka, h, theta, xC, yC, dt)
        zD = self.calculate_zD(m, n, ka, h, theta, xC, yC, dt)
        self.surface[:, 1] = 0.5 * (zT + zD)
        return self.surface

    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot(self.x, self.y)
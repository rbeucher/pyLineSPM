import numpy as np


class River(object):
    
    def __init__(self, surface, nodes, max_node, min_node, precipitation_rate, erodibility, uplift):

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
        self.zT = np.zeros_like(self.surface)
        self.zR = np.zeros_like(self.surface)

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
    
    @staticmethod
    def calculate_D(h, m, n, N, yC, yD):
        hmn = h*m/n
        
        if hmn == 1.0:
            D = -0.5 * N*np.log(yC/yD)
        else:
            D = 0.5 * yD**(1-hmn) * N * 1/(1-hmn) * (1-(yC/yD)**(1-hmn))
        return D
    
    def calculate_zT(self, m=0.5, n=1.0, ka=0.3, h=2.0, dt=0.1, bc=0., xC=0.2):
        """ Calculate and apply erosion in the trunk channel"""
        
        # The Lowest node in the river is assumed to
        # be fixed, the erosion rate is equal to 0.
        # We calculate the value of the erosion rate
        # starting from the second lowest node and moving
        # in the upstream direction.
        
        # We reverse all arrays
        K = self.erodibility[::-1]
        Q = self.get_discharge(ka, h)[::-1]
        zt = self.y[::-1]
        x = self.x[::-1]
        dx = np.diff(x)[0]
        xC = xC * np.max(self.x)
        
        e0 = 1e-6
        edot = e0
        
        # We integrate from the boundary node to the
        # main divide.
        for idx in range(1, len(x)):
            reltol = 1e3
            if x[idx] < xC:
                edot = 0.
            else:
                while reltol > 1e-7:
                    prev = edot
                    edot = -K[idx] * Q[idx]**m * (zt[idx] - zt[idx - 1] + dt * (edot - 0.)) / dx**n
                    reltol = abs((edot - prev) / (edot+1e-8)) # Avoid divide by zero here.
            zt[idx] = zt[idx] + dt * edot
            
        self.zT = zt[::-1]
        return self.zT
    
    def calculate_zR(self, m=0.5, n=1.0, ka=0.3, h=2.0, theta=30, xC=1000, yC=0.1, dt=0.1):
        
        theta = np.radians(theta)
        K = self.erodibility
        P = self.precipitation_rate
        U = self.uplift
        
        xD = np.max(self.x)
        zT = self.zT
        zR = np.zeros_like(self.zT)
        
        hmn = h*m/n

        for ix, x in enumerate(self.x): 
            D = 0.
            yD = 0.5 * ka * x**(h-1.0)
            # Calculate group parameter        
            N = (1.0 / (xD**(h*m) * K[ix] * P[ix]**m * ka**m))**(1/n)
            if yD:
                if hmn == 1.0:
                    D = -0.5 * N*np.log(yC/yD)
                else:
                    D = 0.5 * yD**(1-hmn) * N * 1/(1-hmn) * (1-(yC/yD)**(1-hmn))
            if yD > yC:
                zD = zT[ix] + yC * np.tan(theta) + 2. * D * U[ix]**(1/n)*xD**hmn
            else:
                zD = zT[ix] + yD * np.tan(theta)
        
            edotR = 0.
            if D:
                edotR = K[ix] * P[ix]**m*ka**m * (zD - zT[ix] - yC * np.tan(theta) / (2*D*xD**hmn))**n
            print(edotR)
            zR[ix] = zD + dt * edotR # need to check that...

        return zR
    
    def update_surface(self, m=0.5, n=1.0, ka=0.3, h=2.0, theta=30, xC=1000, yC=100, dt=0.1):
        zT = self.calculate_zT(m, n, ka, h, dt)
        zR = self.calculate_zR(m, n, ka, h, theta, xC, yC, dt)
        self.surface[:, 1] = zT #0.5 * (zT + zR)
        return self.surface

    def plot(self):
        import matplotlib.pyplot as plt
        plt.plot(self.x, self.y)
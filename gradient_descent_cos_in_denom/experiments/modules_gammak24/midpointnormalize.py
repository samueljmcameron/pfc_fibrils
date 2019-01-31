# took this from a stackoverflow discussion on setting the
# midpoint of a density plot or a contour plot to 0. Credit
# goes to Joe Kington. Find the discussion here:
#   https://stackoverflow.com/questions/20144529
#   /shifted-colorbar-matplotlib/20146989#20146989

# this should work for very simple cases, and for imshow,
# contour, or contourf (from matplotlib)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


class MidpointNormalize(Normalize):

    def __init__(self,vmin=None,vmax=None,midpoint=None,clip=False):
        self.midpoint=midpoint
        Normalize.__init__(self,vmin,vmax,clip)

        return

    def __call__(self,value,clip=None):
        # I'm ignoring masked values and all kinds of edge cases to
        # make a simple example...
        x,y = [self.vmin,self.midpoint,self.vmax],[0,0.5,1]
        return np.ma.masked_array(np.interp(value,x,y))
    


if __name__ == "__main__":

    data = np.random.random((10,10))
    
    data = 10*(data-0.8)

    fig,ax = plt.subplots()
    norm = MidpointNormalize(midpoint=0)
    im = ax.imshow(data,norm=norm,cmap=plt.cm.seismic,
                   interpolation='none')
    fig.colorbar(im)
    plt.show()

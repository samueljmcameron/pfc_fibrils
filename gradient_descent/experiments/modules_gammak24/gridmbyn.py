import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec


class GridmByn(object):

    # This class builds the underlying structure of a figure with nrows
    # of subplots per column and ncols of subplots per row, so the total
    # number of subplots is nrows x ncols. Each subplot has three sets
    # of axes, one for a main 3d heatmap plot, one for a slice through
    # the heat map at a constant value of the y coordinate, and one for
    # a slice through the heat map at a constant value of the x
    # coordinate. An example plot (without data) can be produced by just
    # running the exampleplot() method.

    # To access the subplot at gridpoint i,j (which is really a list of
    # the three axes in each subplot), call axarr[i][j]. To access the
    # main axis of a subplot, call axarr[i][j]['main']. To access the
    # constant x (y) slice axis of subplot i,j, call
    # axarr[i][j]['slice_const_x'] (axarr[i][j]['slice_const_y']).

    # The main input of the class is a fig, which can just be generated
    # with plt.figure().

    main_sharex = None
    main_sharey = None
    slice_const_y_sharey = None
    slice_const_x_sharex = None
    
    def __init__(self,fig,nrows=3,ncols=3,wspace=0.1,hspace=0.1):
        self.fig = fig
        self.nrows = nrows
        self.ncols = ncols
        self.wspace = wspace
        self.hspace = hspace
        self.gs0 = self.fig.add_gridspec(nrows,ncols,wspace=self.wspace,
                                         hspace=self.hspace)
        self.axarr = []
        return

    def set_axes_sharing(self,i,j):
        if i == 0 and j == 0:
            self.main_sharex = self.main_sharey = None
            self.slice_const_y_sharey = None
            self.slice_const_x_sharex = None

        else:
            self.main_sharex = self.main_sharey = self.axarr[0][0]['main']
            self.slice_const_y_sharey = self.axarr[0][0]['slice_const_y']
            self.slice_const_x_sharex = self.axarr[0][0]['slice_const_x']

        return


    def make_gssub(self,i,j,nrows=2,ncols=2,width_ratios=[3,1],height_ratios=[1,3],
                   wspace = 0.0,hspace = 0.0):

        gssub = self.gs0[i,j].subgridspec(nrows,ncols,width_ratios=width_ratios,
                                          height_ratios=height_ratios,
                                          wspace=wspace,hspace=hspace)
        
        return gssub

    
    
    def subgrid2x2(self,gssub):

        axes = {}
    
        axes['main'] = self.fig.add_subplot(gssub[1,0],sharex=self.main_sharex,
                                            sharey=self.main_sharey)
        axes['label'] = self.fig.add_subplot(gssub[0,1])
        axes['label'].set_axis_off()
        axes['slice_const_y'] = self.fig.add_subplot(gssub[0,0],sharex=axes['main'],
                                                sharey=self.slice_const_y_sharey)
        axes['slice_const_x'] = self.fig.add_subplot(gssub[1,1],sharey=axes['main'],
                                                sharex=self.slice_const_x_sharex)

    
        axes['slice_const_y'].label_outer()
        axes['slice_const_x'].label_outer()

        return axes

    def subgrid_label(self,gssub,toplabel=r'$\gamma$',botlabel=r'$k_{24}$'):


        
        slabel = self.fig.add_subplot(gssub[0,1])
        slabel.set_axis_off()
        slabel.text(0.1,0.3,toplabel+'\n'+r'$k_{24}=4$')

        return
        

    def build_subplots_subgrid2x2(self,width_ratios=[3,1],height_ratios=[1,3],
                                  wspace=0.0,hspace=0.0):

        for i in range(self.nrows):
            self.axarr.append([])
            for j in range(self.ncols):
                
                gssub = self.make_gssub(i,j,width_ratios=width_ratios,
                                        height_ratios=height_ratios,
                                        wspace=wspace,hspace=hspace)


                self.set_axes_sharing(i,j)
                
                self.axarr[i].append(self.subgrid2x2(gssub))
                if i != 2:
                    plt.setp(self.axarr[i][j]['main'].get_xticklabels(),visible=False)
                    plt.setp(self.axarr[i][j]['slice_const_x'].get_xticklabels(),visible=False)

                if j != 0:
                    plt.setp(self.axarr[i][j]['main'].get_yticklabels(),visible=False)
                    plt.setp(self.axarr[i][j]['slice_const_y'].get_yticklabels(),visible=False)

        return

    def exampleplot(self,width_ratios=[3,1],height_ratios=[1,3],
                    wspace=0.0,hspace=0.0):

        self.build_subplots_subgrid2x2(width_ratios=width_ratios,
                                       height_ratios=height_ratios,wspace=wspace,
                                       hspace=hspace)

        plt.show()

        return
    

# dummy stuff for testing
if __name__ == "__main__":

    fig = plt.figure()
    width = height = 3.487
    fig.set_size_inches(3*width,3*height)


    grid = GridmByn(fig)

    grid.exampleplot()

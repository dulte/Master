import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import time
from math import sqrt

import os
import pickle

from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)


from postcactus.simdir import SimDir
from postcactus import visualize as viz
from postcactus import grid_data as gd

from scipy import ndimage, interpolate

import matplotlib.colors as colors
import seaborn as sns

sns.set_style("darkgrid")
plt.rcParams.update({"font.size": 30})


from ETQuantities import ETQuantities, ETQuantities_gridInterpolator




class ReadQuantities:
    """
    This class wraps ETQuantities and ETQuantities_gridInterpolator and makes it 
    easier to read data from Einstein Toolkit. This class should be used instead of
    using the two classes directly.

    This will only hold one quantity at the time. So when calling read(), the 
    previously held quanity is replaced.

    The user will in most cases only instance this class, then sent it to
    ETInterpolater for use. The only direct use of this class is to make test plots.

    Typical use:

    .. code-block:: python
        g = None
        limits = [250]
        et_q = ReadQuantities([g], it, folder, pickle_folder=pickle_folder, pickle=False, linear=True, limits=limits)
        et_q.test_plot("gxx")
    """
    def __init__(self, geometries, iteration, simulation_folder,  quantity_names=[], pickle=True, pickle_folder="", linear=False, limits=None):
        """

        Parameters
        ----------
        geometry : list of postcactus.grid_data.RegGeom
            List of geometries to be used
            
        iterations: list of ints
            List of iterations to be used

        simulation_folder: str
            Location of the Einstein Toolkit simulation

        quantity_names: list, optional
            List of quantities to read. If nothing is given all quantities will be read (default is [])

        pickle: bool, optional
            Whether or not to pickle the interpolation (default is False)
        
        pickle_folder: str, optional
            Where to save the pickles (default is "")
        
        linear: bool, optional
            Whether to use linear interpolation. If false, a spline is used (default is False)

        limits: list, optional
            Limits where to switch geometries. If nothing is given the the whole geometry is used. (Default is None)


        """
        self.quantity_instances = []
        self.limits = None

        if geometries == None or geometries[0] == None:
            self.geometries = None
            self.quantity_instances.append(ETQuantities(self.geometries, iteration, simulation_folder, quantity_names, pickle, pickle_folder))

            print "[+] Ready to Read Full Grid, with a None geometry"
        else:
            self.geometries = geometries
            self.limits = [abs(g.x1()[0]) for g in self.geometries]
            

            self.geometries = [x for _,x in sorted(zip(self.limits, geometries))]
            for g in self.geometries:
                if linear:
                    self.quantity_instances.append(ETQuantities_gridInterpolator(g, iteration, simulation_folder, quantity_names, pickle, pickle_folder))
                else:
                    
                    self.quantity_instances.append(ETQuantities(g, iteration, simulation_folder, quantity_names, pickle, pickle_folder))
            
            if limits is None or len(limits) != len(self.limits):
                self.limits = sorted(self.limits)
            else:
                self.limits = limits
        
            print "[+] Ready to Read %s Grids, with the Radii %s" %(len(self.limits), self.limits)


    def read(self, name):
        """
        Reads a Einstein Toolkit quantity from HDF5 files for a given iteration
        and geometery from the init function. Takes only the name of the quantity

        Parameters
        ----------
        name : str
            First part of the name of the HDF5 file containg that data
        
        Returns
        -------
        postcactus.grid_data.grid
            Function values of read data, with in the given geometry

        """
        for i, q in enumerate(self.quantity_instances):
            
            if self.geometries == None or len(self.geometries) == 1:
                print "[~] Starting to Read %s on Grid with None geometry" %name
            else:
                print "[~] Starting to Read %s on Grid with Corner %s" %(name,self.geometries[i].x1())

            self.quantity_instances[i].read(name)

            print "[+] Done Reading"
        
    
    def test_plot(self, quantity="alp",binary=False):
        """
        A simple function to test if the module is able to read and interpolate
        a given quantity. The result of the xy plane is then plotted as a
        pcolormesh.

        Parameters
        ----------
        quantity : str
            The quantity the user wants to test plot.

        """

        self.read(quantity)
        """
        try:
            
            self.read(quantity)
        except:
            raise ValueError("[-] Quantity %s not found" %quantity)
        """
        n = 200
        end = 12#20
        q_inter = np.zeros((n,n))
        q_diff = np.zeros((n,n))
        one_dim_plot = np.zeros((n,3))
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        thetas = np.linspace(0,2*np.pi, n)
        n_r = 400
        rs = np.linspace(.5, 250, n_r)

        contours = np.zeros((n, n_r))

    
        #M = 1.0
        M = 0.47656
        
        r_sch = 2*M

        name_tag = quantity
        if quantity == "alp":
            name_tag = r"$\alpha$"
        elif quantity == "gxx":
            name_tag = r"$g_{xx}$"

        for i,r in enumerate(rs):
            for j,theta in enumerate(thetas):
                c = np.array([abs(r*np.cos(theta)), (r*np.sin(theta)),0])
                contours[j,i] = self(c)


        for i in range(n):
            one_dim_plot[i,0] = self([abs(x[i]), 0, 0])
            one_dim_plot[i,1] = self([0,abs(y[i]), 0])
            if quantity == "gxx":
                if binary:
                    r1 = np.abs(3-x[i])
                    r2 = np.abs(-3-x[i])
                    one_dim_plot[i,2] = self([abs(x[i]), 0, 0]) - 0.5*(1+r_sch/(4*abs(r2)))**4 - 0.5*(1+r_sch/(4*abs(r1)))**4
                else:
                    one_dim_plot[i,2] = self([abs(x[i]), 0, 0]) - (1+r_sch/(4*abs(x[i])))**4
                
            elif quantity == "alp":
                one_dim_plot[i,2] = self([abs(x[i]), 0, 0]) - (1-r_sch/abs(4*x[i]))/(1+r_sch/abs(4*x[i]))
                

            

            for j in range(n):
                input = np.array([abs(x[i]), abs(y[j]),0])
                R = sqrt(x[i]**2 + y[j]**2)
                
                q_inter[j,i] = self(input)
                if quantity == "gxx":
                    q_diff[j,i] = q_inter[j,i] - (1+r_sch/(4*R))**4
                elif quantity == "alp":
                    q_diff[j,i] = q_inter[j,i] - (1-r_sch/abs(4*R))/(1+r_sch/abs(4*R))


        plt.pcolormesh(x,y,q_inter)
        plt.colorbar(label=name_tag)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Intesity Map of "+name_tag)
        plt.show()

        plt.contour(x,y,q_inter, 20)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Contour of "+name_tag)
        plt.show()

        plt.plot(x, one_dim_plot[:,0])
        #plt.plot(y, one_dim_plot[:,1])
        if quantity == "gxx":
            plt.title(r"$g_{xx}$")
            plt.ylabel(r"$g_{xx}$")
        elif quantity == "alp":
            plt.title(r"$\alpha$")
            plt.ylabel(r"$\alpha$")
        plt.ylim(-1,10)
        plt.xlabel("r")
        plt.show()


        np.savetxt("%s_x_none.txt" %quantity, x)
        np.savetxt("%s_quantity_none.txt" %quantity, one_dim_plot[:,0])

        """

        plt.plot(x, one_dim_plot[:,0])
        plt.plot(y, one_dim_plot[:,1])
        if quantity == "alp":
            plt.plot(x, (1-1/(2*np.abs(x)))/(1+1/(2*np.abs(x))))
        elif quantity == "gxx":
            plt.plot(x, (1+r_sch/(4*np.abs(x)))**4)
        
        plt.ylim(-1,10)
        plt.show()


        mean_diff = np.zeros(n_r)
        max_diff = np.zeros(n_r)
        for i in range(n_r):
            #norm_value = abs(contours[:,i]-(np.mean(contours[:,i])))/np.mean(contours[:,i])
            if quantity == "gxx":
                norm_value = abs(contours[:,i]-(1+r_sch/(4*abs(rs[i])))**4)
            elif quantity == "alp":
                norm_value = abs(contours[:,i]-(1-r_sch/abs(4*rs[i]))/(1+r_sch/abs(4*rs[i])))
            mean_diff[i] = np.mean(norm_value)
            max_diff[i] = np.max(norm_value)

            plt.plot(thetas,norm_value,label="r=%.2f"%rs[i])
        plt.legend()
        plt.show()


        np.save("r_ET_data", rs)
        np.save("mean_diff_ET_data", mean_diff)
        plt.semilogy(rs, mean_diff)
        plt.title(r"Mean over $\theta$ of $\epsilon(r, \theta) = |g_{xx}^{simulated}-g_{xx}^{analytical}|$")
        plt.xlabel("r")  
        plt.ylabel(r"$mean_\theta (\epsilon(r, \theta))$")  
        plt.show()
        
        plt.semilogy(rs, max_diff)    
        plt.title(r"Max over $\theta$ of $\epsilon(r, \theta) = |g_{xx}^{simulated}-g_{xx}^{analytical}|$")
        plt.xlabel("r")  
        plt.ylabel(r"$max_\theta (\epsilon(r, \theta))$")
        plt.show()

        plt.pcolormesh(x,y,np.abs(q_diff), norm=colors.LogNorm())
        plt.colorbar(label=name_tag)
        if quantity == "gxx":
            plt.title(r"$|g_{xx}^{simulated}-g_{xx}^{analytical}|$")
        elif quantity == "alp":
            plt.title(r"$|\alpha^{simulated}-\alpha^{analytical}|$")
        plt.show()

        plt.contour(x,y,np.abs(q_diff),norm=colors.LogNorm())
        if quantity == "gxx":
            plt.title(r"$|g_{xx}^{simulated}-g_{xx}^{analytical}|$")
        elif quantity == "alp":
            plt.title(r"$|\alpha^{simulated}-\alpha^{analytical}|$")
        plt.show()
        
        """

        if quantity == "gxx" or quantity == "alp":
            """
            plt.plot(x, one_dim_plot[:,2])
            plt.ylim(-1,0.5)
            plt.show()
            """

            plt.semilogy(x[1:-1], np.abs(one_dim_plot[1:-1,2]))
            plt.xlabel("radius")
            plt.ylabel("|Simulated-Analytical|")
            plt.title("|Simulated-Analytical| for "+name_tag)
            print np.mean(np.abs(one_dim_plot[1:-1,2]))
            print np.max(np.abs(one_dim_plot[1:-1,2]))


            
            print "Mean: " + str(np.mean( np.abs( one_dim_plot[np.abs(x)>0.5,2])[1:-1] ))
            print "Max: " + str(np.max( np.abs( one_dim_plot[np.abs(x)>0.5,2])[1:-1] ))

            plt.show()

    """
    def _bound(self, val, axis):
        if axis != "y":
            return abs(val)
        else:
            return val
    """
            


    def __call__(self, coords):
        """
        Returns the the value of the read quantity at the given coordinate.
        Will take case of which geometry to read from.

        Parameters
        ----------
        coords: list
            An 3D array/list of coodinates
        
        Returns
        -------
        float
            The funciton value at the coordinates.


        """

        

        if self.geometries == None:
            return self.quantity_instances[0](coords)
        else:
            for i,q in enumerate(self.quantity_instances):
                limit = self.limits[i]
                if coords[0] > limit or coords[1] > limit or coords[2] > limit:
                    continue
                else:  
                    return self.quantity_instances[i](coords)
            
            return self.quantity_instances[-1](coords)


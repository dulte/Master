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





class ETQuantities(object):
    """
    Class for reading Einstein Toolkit data and interpolating it. This class will use
    a spline interpolator. Use the wrapper ReadQuantities instead of using this
    class directly.
    """
    def __init__(self, geometry, iteration, simulation_folder, quantity_names=[], pickle=True, pickle_folder=""):
        """

        Parameters
        ----------
        geometry : postcactus.grid_data.RegGeom
            Geometry to be used
            
        iterations: int
            Iteration to be used

        simulation_folder: str
            Location of the Einstein Toolkit simulation

        quantity_names: list, optional
            List of quantities to read. If nothing is given all quantities will be read (default is [])

        pickle: bool, optional
            Whether or not to pickle the interpolation (default is False)
        
        pickle_folder: str, optional
            Where to save the pickles (default is "")


        """
        if quantity_names == []:
            self.quantity_names = ["alp"]#, "betax", "betay", "betaz",
                    #"gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
                    #"kxx", "kxy", "kxz", "kyy", "kyz", "kzz"]
        else:
            self.quantity_names = quantity_names

        self.quantities = {}
        self.sd = SimDir(simulation_folder)
        self.geo = geometry
        self.iteration = iteration
        self.pickle = pickle

        self.loaded_spline = None
        if not(geometry is None):
            corner = self.geo.x1()
            self.folder = pickle_folder
            self.filename = "et_quantities_%d_%d_%d" %(corner[0], corner[1], corner[2])
        else:
            self.folder = ""
            self.filename = ""

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
        start_time = time.time()
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)
        self.name = name
        

        if self._check_pickle(name):
            print "[~] Pickle Found. Reading Pickle."
            self.load_pickle(name)

        else:
            print "[~] Pickle Not Found. Reading From ET Files."

            grid = self._read_quantity(name, self.geo, self.iteration)
            if self.geo is None:
                self.loaded_spline = grid
                #self.loaded_spline = self.interpolate_none_geo(grid)
            else:
                self.loaded_spline = grid.spline(order=3, mode="nearest")

        if self.pickle == True and not os.path.exists(full_path) and not(self.geo is None):
            self.pickle_quantities(name)

        print "Read All Quantities in %s seconds" %(time.time()-start_time)
        return self.quantities





    def _check_pickle(self, name):
        """
        Checks if a pickle file exists.

        Parameters
        ----------
        name: str
            Name of the pickle file
        
        Returns
        -------
        bool
            Whether the pickle file exists or not


        """
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)
        if not pickle:
            return False
        elif os.path.exists(full_path):
            return True
        else:
            return False


    def _read_quantity(self,quantity, geometry, iteration, dimentions=3, order=4):
        """
        Reads a Einstein Toolkit quantity from HDF5 files for a given iteration
        and geometery. If the dimensions is 2 the xy plane is returned, if
        3 xyz is return, else a error is raised.

        Parameters
        ----------
        quantity : str
            First part of the name of the HDF5 file containg that data
        geometry : postcactus.grid_data.RegGeom
            The geometry at which the data should be read and interpolated
        iteration : int
            The timestep at which the data should be read
        dimentions : int, optional
            The dimension of which the data should be read (default is 3)
        order : int, optional
            The order of the interpolation (default is 4)

        Raises
        ------

        ValueError
            If the dimensions are not 2 or 3
        ValueError
            If the quantity could not be read/wrong name of quantity


        Returns
        -------
        postcactus.grid_data.grid
            Function values of read data, with in the given geometry


        """


        if dimentions == 2:
            try:
                grid = self.sd.grid.xy.read(quantity, iteration, geom=geometry, order=order)
            except:
                raise ValueError("Quantity %s could not be read!" %quantity)

        elif dimentions == 3:
            
            


            
            try:
                
                grid = self.sd.grid.xyz.read(quantity, iteration, geom=geometry, order=order)
                
            except:
                raise ValueError("Quantity %s could not be read!" %quantity)

            
            
        else:
            raise ValueError("Number of dimentions should be 2 or 3!")

        print "[+] %s Successfully Read" %quantity
        return grid

    




    def pickle_quantities(self, name):
        """
        Makes a pickle of the read interpolated file.

        Parameters
        ----------
        name: str
            The name of the pickle file
        
    
        """
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)
        if self.pickle:
            f = open(full_path, "w")
            pickle.dump(self.loaded_spline, f)

    def load_pickle(self, name):
        """
        Loads a pickles file
        Parameters
        ----------
        name: str
            Name of file
        

        """
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)

        f = open(full_path, "r")
        self.loaded_spline = pickle.load(f)


    def test_plot(self, quantity="alp"):
        """
        A simple function to test if the module is able to read and interpolate
        a given quantity. The result of the xy plane is then plotted as a
        pcolormesh.

        Parameters
        ----------
        quantity : str
            The quantity the user wants to test plot.

        """

        
        try:
            q = self.read(quantity)
        except:
            raise ValueError("[-] Quantity %s not found" %quantity)

        n = 100
        end = 50
        q_inter = np.zeros((n,n))
        one_dim_plot = np.zeros(n)
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        for i in range(n):
            one_dim_plot[i] = self([abs(x[i]), 0, 0])

            for j in range(n):
                input = np.array([abs(x[i]), abs(y[j]),0])

                q_inter[j,i] = self(input)


        plt.pcolormesh(x,y,q_inter)
        plt.colorbar()
        plt.show()

        plt.contour(x,y,q_inter, 20)
        plt.show()

        plt.plot(x, one_dim_plot)
        plt.ylim(-1,10)
        plt.show()


    def __call__(self, coords, output=None):
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

        q = self.loaded_spline
        if q is None:
            raise ValueError("[-] No Quantity Loaded")

        if self.name in ["alp", "gxx", "gyy", "gzz"]:
            if sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) > 280:
                return 1
            #elif sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) < 0.0005:
            #    return q([0.0005,0,0])


        if self.geo == None:
            return q(coords)

        if (len(coords) != len(q.data.shape)):
            raise ValueError('Dimension mismatch with sampling coordinates.')
        #

        ind     = [[(c - c0)/dx] for c,c0,dx in zip(coords, q.x0, q.dx)]
        prefilt = (q.order < 2)
        res = ndimage.map_coordinates(q.data, ind, prefilter=prefilt,
             mode=q.mode, cval=q.outside_value, order=q.order,
             output=output)
        if output is None:
            return res
        return output




class ETQuantities_gridInterpolator(ETQuantities):
    """
    Class for reading Einstein Toolkit data and interpolating it. This class will use
    a linear interpolator. Use the wrapper ReadQuantities instead of using this
    class directly.
    """
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
        start_time = time.time()
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)
        self.name = name

        if self._check_pickle(name):
            print "[~] Pickle Found. Reading Pickle."
            self.load_pickle(name)

        else:
            print "[~] Pickle Not Found. Reading From ET Files."

            grid = self._read_quantity(name, self.geo, self.iteration)
            self.loaded_spline = self._get_interpolator(grid)

        if self.pickle == True and not os.path.exists(full_path):
            self.pickle_quantities(name)

        print "Read All Quantities in %s seconds" %(time.time()-start_time)
        return self.quantities

    def _get_interpolator(self, grid):
        """
        Creates a linear interpolation given a read grid of data

        Parameters
        ----------
        grid: postcactus.grid_data.grid
            Grid of Einstein Toolkit data
        
        Returns
        -------
        Objects
            The callable interpolation object


        """

        #shape = self.geo.shape()
        
        shape = grid.data.shape
        min_corner = self.geo.x0()
        max_corner = self.geo.x1()
        
        x = np.linspace(min_corner[0], max_corner[0], shape[0])
        y = np.linspace(min_corner[1], max_corner[1], shape[1])
        z = np.linspace(min_corner[2], max_corner[2], shape[2])

        return interpolate.RegularGridInterpolator((x,y,z), grid.data, bounds_error=False, fill_value=None)

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
        
        #if self.name in ["alp", "gxx", "gyy", "gzz"]:
            #if sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) > 800:
            #    return 1
            #elif sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) < 1.5:
            #    return 0
        
        return self.loaded_spline(coords)
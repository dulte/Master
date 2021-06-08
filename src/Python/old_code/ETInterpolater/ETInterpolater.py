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



"""
TODO:
- Add change of origin of coord system [~]
- Guard against inf and nan in x,y,z when getting the values [+]
- Conversion between LORENE and ET units [+]
- Take care of symmetries from ET [~]
    - As of now this only takes abs, but it needs to get symmetry axis --if any.
- Make sure boudaries are handled correctly! [-]
    - PostCactus cannot extrapolate outside the ET grid.
      So this grid must be larger than highest coll point

- Is the number of coll points hard coded?
"""


class ReadQuantities:
    def __init__(self, geometries, iteration, simulation_folder,  quantity_names=[], pickle=True, pickle_folder="", linear=False, limits=None):
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


    def _bound(self, val, axis):
        if axis != "y":
            return abs(val)
        else:
            return val

            


    def __call__(self, coords):

        

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









class ETQuantities(object):
    def __init__(self, geometry, iteration, simulation_folder, quantity_names=[], pickle=True, pickle_folder=""):
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

    def interpolate_none_geo(self, grid):
        """
        points = []
        values = []
        grid_coords = grid.coords()

        for i,compdata in enumerate(grid_coords):
            for j,d in enumerate(compdata):
                coords1d = d.coords1d()
                l = len(coords1d[0])
                
                points += [[coords1d[0][k],coords1d[1][k],coords1d[2][k]] for k in range(l)]
        for point in points:
            values.append(grid(point))

        
        min_corner = 0
        max_corner = 250
        shape = 250
        points = np.array(points)
        values = np.array(values)
        
        points = np.where(np.abs(points)< 1e-10, 1e-10, points)
        values = np.where(np.abs(values)< 1e-10, 1e-10, values)
        
        
        #x = np.linspace(min_corner, max_corner, shape)
        #y = np.linspace(min_corner, max_corner, shape)
        #z = np.linspace(min_corner, max_corner, shape)

        x = points[:,0]
        y = points[:,1]
        z = np.array(points[:,2]) + np.random.rand(len(x))*1e-8
        """

        # put the available x,y,z data as a numpy array
        points = np.array([[ 27.827,  18.53 , -30.417], [ 24.002,  17.759, -24.782],[ 22.145,  13.687, -33.282], [ 17.627,  18.224, -25.197],[ 29.018,  18.841, -38.761], [ 24.834,  20.538, -33.012],[ 26.232,  22.327, -27.735], [ 23.017,  23.037, -29.23 ],[ 28.761,  21.565, -31.586], [ 26.263,  23.686, -32.766]])
        # and put the moisture corresponding data values in a separate array:
        values = np.array([0.205,  0.197,  0.204,  0.197,  0.212, 0.208,  0.204,  0.205, 0.211,  0.215])
        request = np.array([[25, 0, -30], [27, 20, -32]])
        print points.shape
        print values.shape
        data = interpolate.LinearNDInterpolator(points, values)
        
        #print(data(np.array([[10,11,10], [-1, 0, 0]])))
        print(data(request))
        exit()

        return data#interpolate.RegularGridInterpolator((x,y,z), data, bounds_error=False, fill_value=None)



    def pickle_quantities(self, name):
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)
        if self.pickle:
            f = open(full_path, "w")
            pickle.dump(self.loaded_spline, f)

    def load_pickle(self, name):
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
    def read(self, name):
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
        best results is 200 and 1.5
        """
        
        if self.name in ["alp", "gxx", "gyy", "gzz"]:
            if sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) > 800:
                return 1
            #elif sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2) < 1.5:
            #    return 0
        
        return self.loaded_spline(coords)




class ETInterpolater:
    """
    A class used to read the 3+1 quantities simulated with Einstein Toolkit,
    interpolate them, then use C code build around LORENE to find the collocation
    points needed to do a spectral transformation. It can then find the function
    vaules at these points, and call the C code to do the spectral transformation.

    Most user will only need to use make_geometry/make_positive_geometry and analyse_bbh.
    A typical run example will be

    .. code-block:: python

        folder = "/some_location/simulations/bbh_3D"
        inter = ETInterpolater(folder, 2)

        g = inter.make_positive_geometry([100, 100, 100], 200)
        it = [0, 128, 256] #And more

        inter.analyse_bbh(g, quantity="", it, test=True)
    """

    def __init__(self, dir, nb_bodies=2):
        """
        Parameters
        ----------

        dir : str
            The directory of the Einstein Toolkit results
        nb_bodies : int, optional
            Number of objects in the Einstein Toolkit simulation (default set to 2)

        """
        self.sd = SimDir(dir)
        self.nb_bodies = nb_bodies

        self.dir = dir

        self.xlim = [0,0]
        self.ylim = [0,0]
        self.zlim = [0,0]


    def read_ET_quantity(self,quantity, geometry, iteration, dimentions=3, order=4):
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


    def make_geometry(self, corner, n_pts):
        """
        Make a geometry needed to read Einstein Toolkit quantities.
        Given the bottom left corner, this will return a geometry using
        a second corner which is a mirror of the given corner, mirrored
        around the origin.

        Parameters
        ----------
        corner : list
            Bottom left corner of the geomtery
        n_pts : int
            Number of points in the geomtery


        Raises
        ------
        ValueErrror
            If the corner is not 2 or 3 dimensional

        Returns
        -------
        postcactus.grid_data.RegGeom
            The geometry


        """
        if len(corner) > 3 or len(corner) < 2:
            raise ValueError("Wrong Number of Corners! User 2 or 3!")

        corner1 = [-i for i in corner]
        self.xlim = [corner1[0], corner[0]]
        self.ylim = [corner1[1], corner[1]]

        if len(corner) == 3:
            self.zlim = [corner1[2], corner[2]]

        geo = gd.RegGeom([n_pts]*len(corner), corner, x1=corner1)
        print "[+] Geomentry Successfully Made"
        return geo


    def make_semipositive_geometry(self, corner, n_pts):
        if len(corner) > 3 or len(corner) < 2:
            raise ValueError("Wrong Number of Corners! User 2 or 3!")

        corner1 = [abs(i) for i in corner]
        self.xlim = [-abs(corner1[0]), abs(corner[0])]
        self.ylim = [-abs(corner1[1]), abs(corner[1])]

        start_corner = [0, -abs(corner1[1])]

        if len(corner) == 3:
            self.zlim = [abs(corner[2]), abs(corner[2])]
            start_corner.append(0)

        geo = gd.RegGeom([n_pts]*len(corner), start_corner, x1=corner1)
        print "[+] Geomentry Successfully Made"
        return geo




    def make_positive_geometry(self, corner, n_pts):
        """
        Make a geometry needed to read Einstein Toolkit quantities.
        This will make a geometry with one corner at the origin, and
        the other at the given corner. This is for cases where the
        Einstein Toolkit simulation is done with symmetry in xyz!

        Parameters
        ----------
        corner : list
            Bottom left corner of the geomtery
        n_pts : int
            Number of points in the geomtery


        Raises
        ------
        ValueErrror
            If the corner is not 2 or 3 dimensional

        Returns
        -------
        postcactus.grid_data.RegGeom
            The geometry


        """

        if len(corner) > 3 or len(corner) < 2:
            raise ValueError("Wrong Number of Corners! User 2 or 3!")

        corner1 = [abs(i) for i in corner]
        self.xlim = [-abs(corner[0]), abs(corner[0])]
        self.ylim = [-abs(corner[1]), abs(corner[1])]

        if len(corner) == 3:
            self.zlim = [abs(corner[2]), abs(corner[2])]

        geo = gd.RegGeom([n_pts]*len(corner), [0]*len(corner), x1=corner1)
        print "[+] Positive Geomentry Successfully Made"
        return geo



    def make_test_plot(self, quantity):
        """
        A simple function to test if the module is able to read and interpolate
        a given quantity. The result of the xy plane is then plotted as a
        pcolormesh.

        Parameters
        ----------
        quantity : str
            The quantity the user wants to test plot.

        """
        g = self.make_geometry([-50,-50], 400)
        q = self.read_ET_quantity(quantity, g, 0, dimentions=2)

        n = 300
        end = 10
        q_inter = np.zeros((n,n))
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        for i in range(n):
            for j in range(n):
                q_inter[j,i] = q(np.array([self.desymmetrize_coord(x[i]), self.desymmetrize_coord(y[j])]))


        plt.pcolormesh(x,y,q_inter)
        plt.colorbar()
        plt.show()


    def make_test_bbh_plot(self, quantity, p1, p2, r1, r2):
        """
        A simple function to test if the module is able to read and interpolate
        a given quantity with two bodies. The code wil apply the splitt function.
        The result of the xy plane is then plotted as a pcolormesh and a contour plot.

        Parameters
        ----------
        quantity : str
            The quantity the user wants to test plot.

        """
        #g = self.make_geometry([-50,-50], 400)
        g = self.make_positive_geometry([50,50], 400)
        q = self.read_ET_quantity(quantity, g, 0, dimentions=2)

        n = 600
        end = 10
        q_inter = np.zeros((n,n))
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        th = np.linspace(0,2*np.pi,100)

        scaling_factor = 4.0
        bbh_distance = np.linalg.norm(p1-p2)

        for i in range(n):
            for j in range(n):
                q_inter[j,i] = q(np.array([self.desymmetrize_coord(x[i]), self.desymmetrize_coord(y[j])]))
                q_inter[j,i] *= self.split_function(x[i], y[j], 0, p1, p2, bbh_distance/scaling_factor, bbh_distance/scaling_factor)


        xx = p1[0] + r1*np.cos(th)
        yy = p1[1] + r1*np.sin(th)
        plt.plot(xx, yy, linestyle=':', color='black')

        xx = p2[0] + r2*np.cos(th)
        yy = p2[1] + r2*np.sin(th)
        plt.plot(xx, yy, linestyle=':', color='black')

        plt.plot(p1[0], p1[1], "r*")
        plt.plot(p2[0], p2[1], "r*")

        plt.pcolormesh(x,y,q_inter)
        plt.colorbar()
        plt.show()

        plt.contour(x,y,q_inter)
        plt.show()


    def get_values_at_coll_points(self, interpolated_quantity, c_path="../../C/", smooth=True, bh_pos=[0,0,0], bh_rad=0,bh_pos2=[0,0,0], bh_rad2=0, scaling_factor=4.0, test=False):
        """
        One of the main functions of the module. This function takes in a interpolation
        function of a quantity. It will then use LORENE to find the collocation points.
        It will then go though all the collocation points, find the fuction value
        at that point, then use the splitting function on the value (given smooth is true)

        Parameters
        ----------
        interpolated_quantity : postcactus.grid_data.grid
            Interpolation function for the given quantity
        smooth : bool, optional
            Whether or not the splitting function will be applied (default is True)
        bh_pos : list, optional
            The position of the first BH (default is [0,0,0])
        bh_rad : float, optional
            The radius of the first BH (dafault is 0)
        bh_pos2 : list, optional
            The position of the second BH (default is [0,0,0])
        bh_rad2 : float, optional
            The radius of the second BH (dafault is 0)
        scaling_factor : float, optinal
            The scaling factor used in the splitting funciton (default is 4.0)
        test : bool, optional
            If True this will make a test plot

        Returns
        -------
        Dict :
            The values as a dict with the same structure as the collocation points
        List :
            A flat list of the calculated values
quantities[quantity]


        """
        q = interpolated_quantity

        xx, yy, zz = self.get_coll_points(origin=bh_pos, c_path=c_path)

        bbh_distance = np.linalg.norm(np.array(bh_pos)-np.array(bh_pos2))


        r_test = []
        x_test = []
        y_test = []
        z_test = []
        q_test = []

        print "[+] Coll Points Successfully Formatted"

        values = xx[:]

        flatten_values = []


        if test:
            return values, self.flatten_dict(values)

        for index,domain in enumerate(xx):
            for k in range(len(domain.keys())):
                for j in range(len(domain[k].keys())):
                    for r in range(len(domain[k][j])):
                        x = xx[index][k][j][r]
                        y = yy[index][k][j][r]
                        z = zz[index][k][j][r]


                        
                        inter_q = q([self.bound_coord(x, "x"),self.bound_coord(y, "y"),self.bound_coord(z, "z")])

                        if smooth:
                            #inter_q *= self.super_gaussian(x,y,z, bh_pos, bh_rad)
                            #inter_q *= self.super_gaussian(x,y,z, bh_pos2, bh_rad2)
                            inter_q *= self.split_function(x, y, z, bh_pos, bh_pos2, bbh_distance/scaling_factor, bbh_distance/scaling_factor)

                        values[index][k][j][r] = inter_q

                        if inter_q != inter_q:
                            print index, k, j, r, inter_q
                        if test:

                            r_test.append(np.sqrt(x**2 + y**2 + z**2))
                            x_test.append(x)
                            y_test.append(y)
                            z_test.append(z)
                            q_test.append(inter_q)
                            #if inter_q == 0:
                            #    print x,y,z

                        flatten_values.append(inter_q)

        print "[+] Interpolated Quantety Successfully Found"

        if test:
            plt.plot(r_test, q_test, ".")
            plt.show()

        return values, flatten_values


    def bound_coord(self, coord, tp):
        """
        Bounds the coordinates, so that other functions don't try to get
        function values outside of defined area. It also handles infinite
        bounderies which is sometime needed. It will also desymmeterize
        the coordinate (see the desymmeterize_coord function)

        Parameters
        ----------
        coord : float
            The value of the coordinate
        tp : str
            Which axis: x, y or z

        Raises
        ------
        ValueError
            If tp is not x, y or z

        Returns
        -------
        float
            The bounded and desymmeterized coodinate value

        Notes
        -----
        The infinite bounds have to be handled better. It will
        now use the outer coordinate - 5, which is garantied
        to lead to discontinuities...
        """


        if coord == np.inf or coord != coord:
            coord = 50000
            """
            if tp == "x":
                coord = self.xlim[0] - 5
            elif tp == "y":
                coord = self.ylim[0] - 5
            elif tp == "z":
                coord = self.zlim[0] - 5
            else:
                raise ValueError("%s is not an axis" %tp)
            """
        elif coord == -np.inf:# or coord == -np.nan:
            coord = -50000
            """
            if tp == "x":
                coord = self.xlim[1] + 5
            elif tp == "y":
                coord = self.ylim[1] + 5
            elif tp == "z":
                coord = self.zlim[1] + 5
            else:
                raise ValueError("%s is not an axis" %tp)
            """

        return self.desymmetrize_coord(coord)
        """
        if tp != "y":
            return self.desymmetrize_coord(coord)
        else:
            return coord
        """
    def desymmetrize_coord(self, coord):
        """
        Function for taking care of the symmetries used in Einstein Toolkit.
        This is a preliminary function, that assumes that Einstein Toolkit
        uses symmetry in both x, y and z.

        Parameters
        ----------
        coord : float
            The value of the coordinate (point) that needs to be desymmeterized

        Returns
        -------
        float
            Desymmeterized coordinate point.

        Notes
        -----
        This is a temperary function. It assumes total symmetery, and therefore
        only uses a abs() to desymmeterize. Parameters to choose symmeterization
        axis will be added later.
        """
        return abs(coord)



    def flatten_dict(self, d):
        """
        Takes a dict formatted as the output from LORENE and flattens it.

        Parameters
        ----------
        d : dict
            A dict with the values at the collocation points. The dict
            is formatted as the output from LORENE.

        Returns
        -------
        list
            A flatted list with the values from the dict

        Notes
        -----
        To get back the value from the flatted list one can use
        ``` d[l][k][j][i] =  (sum_{0 <= m < l} nr[m]*nt[m]*np[m]) + k*nt[l]*nr[l] + j*nr[l] + i```

        """
        flatten_values = []
        for index,domain in enumerate(d):
            for k in domain.keys():
                for j in domain[k].keys():
                    for r in range(len(domain[k][j])):
                        flatten_values.append(d[index][k][j][r])

        return flatten_values



    def write_values_to_file(self, values, file):
        """
        Writes dict to file

        Parameters
        ----------
        values : dict
            Dict with the same formatting as the output of LORENE
        file : str
            Filename of the file the user wants to save to

        Notes
        -----
        Outdated and should not be used anymore

        """
        with open(file, "w") as f:
            for index,domain in enumerate(values):
                f.write("domain %d: \n \n" %index)
                for k in domain.keys():
                    f.write("k = %d \n" %k)
                    for j in domain[k].keys():
                        f.write("j = %d " %j)
                        for r in range(len(domain[k][j])):
                            f.write("%s " %domain[k][j][r])
                        f.write("\n")
                    f.write("\n")

    def write_flatten_values_to_file(self, values, it, body, file):
        """
        Writes flatten list to file. Adds which body the user is saving,
        the total number of bodies and the iteration/timestep to the
        list before saving

        Parameters
        ----------
        values : list
            Flatten list with the same formatting as the output of LORENE
        it : int
            The iteration/timestep
        body : int
            Which body this is.
        file : str
            Filename of the file the user wants to save to

        """
        values = [body, self.nb_bodies, it] + values
        with open(file, "w+") as f:
            for i in values:
                print_nb = i
                if print_nb != print_nb:
                    print_nb = np.inf
                f.write("%.10f " %print_nb)

        print "[+] File Successfully Written"
        #np.savetxt(file, np.array(values))


    def LORENE_read(self, filename, c_path="../../C/", origin=[0,0,0], body=1, it=0):
        """
        Function responsible of communicating with the C code, which will read
        the flatten array, make the spectral transformation and save the
        Gyoto files.

        Parameters
        ----------
        filename : str
            Name of file users used to save flatten list (not used any longer)
        origin : list, optional
            Origin of the spectral grid LORENE will make (default is [0,0,0])
        body : int, optional
            Which body this is (default is 1)
        it : int, optional
            The iteration/timestep (default is 0)

        Raises
        ------
        IOError
            If the function fails to contact the C code, or if the C code crashes.

        Notes
        -----
        The C code is changed to uses generic insted of a user given filename.
        This means that the filename is irrelevant. This parameter will be
        removed.

        """
        p = subprocess.Popen("%s/get_points %s %s %s 1 %s %s" %(c_path,origin[0], origin[1], origin[2], body, it), stdout=subprocess.PIPE, shell=True)
        #p = subprocess.Popen("./get_points %s %s %s 1" %(origin[0], origin[1], origin[2]), stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        if p_status != 0:
            raise IOError("Could not read LORENE C code!")
        else:
            print "[+] LORENE Successfully Read and Saved tquantities[quantity]he File"



    def get_coll_points(self,c_path="../../C/", origin=[0,0,0], body=1, it=0):
        """
        Function responsible of communicating with the C code to get the
        collocation points, then make the output of the C code to three
        dicts (one for each axis).

        The structure of the dict (and LORENE output) is:

        - ``dict[l]`` is the l'th domain
        - ``dict[l][k]`` is the k'th phi point(s)
        - ``dict[l][k][j]`` is the k'th theta point(s)
        - ``dict[l][k][j][i]`` is the i'th r point

        Parameters
        ----------
        origin : list, optional
            Origin of the spectral grid LORENE will make (default is [0,0,0])
        body : int, optional
            Which body this is (default is 1)
        it : int, optional
            The iteration/timestep (default is 0)

        Returns
        -------
        list, list, list
            Returns a list of dicts for x, y and z with the collocation points.

        Notes
        -----
        This will call all the cleanup function, so this is the only function
        the user need to use.

        """
        p = subprocess.Popen("%s/get_points %s %s %s 0 %s %s" %(c_path,origin[0], origin[1], origin[2], body, it), stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        if p_status != 0:
            raise IOError("Could not read LORENE C code! %s" %err)
        else:
            print "[+] LORENE Successfully Gave the Coll Points"

        with open("printout.txt","w") as f:
            f.write(output)
        #print output
        x_dicts, y_dicts, z_dicts = self.clean_coll_points_xyz(output)
        return x_dicts, y_dicts, z_dicts



    def clean_coll_points_xyz(self,s):
        """
        Takes a string containing the output of the LORENE code and returns
        three dicts (x,y,z) giving the coord position of each coll point.


        Parameters
        ----------
        s : str
            The string given by C code, containg the whole output as one string

        Returns
        -------
        list, list, list
            Returns a list of dicts for x, y and z with the collocation points.

        Notes
        -----
        Very bad practis is used in this function, so might need to be rewritten.
        """
        corrds = []
        index = 0
        for index, sub_s in enumerate(s):
                if sub_s == "+" and s[index+1] == "\n":
                    break
        s = s[index+1:]
        for index, sub_s in enumerate(s):
                if sub_s == "+" and s[index+1] == "\n":
                    break


        x_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        for index, sub_s in enumerate(s):
                if sub_s == "+" and s[index+1] == "\n":
                    break


        y_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        for index, sub_s in enumerate(s):
                if sub_s == "+" and s[index+1] == "\n":
                    break

        z_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        return x_dicts, y_dicts, z_dicts


    def clean_coll_points(self,s):
        """
        Takes the string containing the coll points for one coord (x,y or z),
        and returns a list with a dict for each domain.

        Parameters
        ----------
        s : str
            The string containing the collocation points for one coord (x,y or z)

        Returns
        -------
        list
            Returns a list of dicts for the one coordinate with the collocation points

        """
        domains = []
        s = s.split("\n")
        while len(s) > 0:
            for index, sub_s in enumerate(s):

                if "*** Tbl 3D" in sub_s:

                    break
            domains.append(s[:index-1])
            s = s[index+1:]

        domains = domains[1:]
        domain_dicts = []
        for domain in domains:
            domain_dicts.append(self.get_domain_dict_from_string(domain))

        return domain_dicts

    def get_domain_dict_from_string(self,domain_string):
        """
        Takes a string for one domain of one coord
        and returns the dict for that domain.

        Parameters
        ----------
        domain_string : str
            A string containg the collocation points for one domain and one coord,
            formatted as the output from LORENE

        Returns
        -------
        dict
            Dict with the collocation points for that domain.
        """
        domain_dict = {}
        current_k = None
        for sub_s in domain_string:

            if "k =" in sub_s:
                current_k = int(sub_s.split()[2])
                domain_dict[current_k] = {}
            elif "j = " in sub_s:
                current_j = int(sub_s.split()[2])

                domain_dict[current_k][current_j] = [float(j) for j in sub_s.split()[4:]]

        return domain_dict


    def super_gaussian(self, x, y, z, bh_pos, rad, I=1, n=10):
        """
        Function used to suppress regions of the function grid.

        Parameters
        ----------quantities[quantity]
        x : float
            The x coordinate
        y : float
            The y coordinate
        z : float
            The z coordinate
        bh_pos : list
            The position of the BH the user wants to suppress
        rad : float
            The radius of the suppression
        I : float, optional
            The intesity of the suppression, so the defaul of the function far away
            from the BH (default is 1)
        n : int
            How fast the suppression is (default is 10)

        Returns
        -------
        float
            The suppression factor
        """
        r = np.sqrt((bh_pos[0]-x)**2 + (bh_pos[1]-y)**2 + (bh_pos[2]-z)**2)/float(rad)
        return 1-I*np.exp(-2*r**n)


    def split_function(self, x, y, z, posBH1, posBH2, R1, R2):
        r"""
        C^2 function `f` that takes the value 1 in the vicinity of BH1,
        0 in the vicinity of BH2 and 1/2 far from BH1 and BH2

        INPUT:

        - ``x``, ``y``, ``z`` -- Cartesian coordinates `(x,y,z)` of the point
          where the function `f` is to be evaluated
        - ``posBH1`` -- 3-tuple of Cartesian coordinates specifying the
          location of BH1
        - ``posBH2`` -- 3-tuple of Cartesian coordinates specifying the
          location of BH2
        - ``R1`` -- radius of ball around BH1 where `f(x,y,z) = 1`
        - ``R2`` -- radius of ball around BH2 where `f(x,y,z) = 0`

        OUTPUT:

        - value of `f(x,y,z)`

        """
        x1, y1, z1 = posBH1
        r1 = np.sqrt((x - x1)**2 + (y - y1)**2 + (z - z1)**2)
        if r1 < R1:
            return 1
        x2, y2, z2 = posBH2
        r2 = np.sqrt((x - x2)**2 + (y - y2)**2 + (z - z2)**2)
        if r2 < R2:
            return 0
        D = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        a = D / (R1 + R2)
        A1 = a * R1
        if r1 < A1:
            xx = (r1 - R1) / (A1 - R1)
            return -3*xx**5 + 7.5*xx**4 - 5*xx**3 + 1  # S2_1
        A2 = a * R2
        if r2 < A2:
            xx = (r2 - R2) / (A2 - R2)
            return 3*xx**5 - 7.5*xx**4 + 5*xx**3  # S2_2
        return 0.5



    def read_bbh_diag(self):
        """
        Reads the iterations saved, the positions and radii of the black holes

        Returns
        -------
        np.array, np.array([[list, list, list], [list, list, list]]), np.array([list, list])
            First all the iterations save are returns. Secondly an array with the posision
            of the two BHs over time. Lastly the an array with the radii of the BHs
            over time.
        """
        itbh1, tbh1, xbh1, ybh1, zbh1, rbh1 = np.loadtxt("%s/output-0000/bbh3d/BH_diagnostics.ah1.gp" %self.dir,
                                                         usecols=(0,1,2,3,4,7),
                                                         unpack=True)

        itbh2, tbh2, xbh2, ybh2, zbh2, rbh2 = np.loadtxt("%s/output-0000/bbh3d/BH_diagnostics.ah2.gp" %self.dir,
                                                         usecols=(0,1,2,3,4,7),
                                                         unpack=True)


        return np.array(itbh1), np.array([[xbh1, ybh1, zbh1], [xbh2, ybh2, zbh2]]), np.array([rbh1, rbh2])





    def analyse_bbh(self, geometry, ETquantities, iterations, c_path="../../C/", result_path="./", do_gyoto_converstion=True, quantities=[], test=False, split=True, scaling_factor=4.0):
        """
        The main function of the code. This will read and interpolate all quantities,
        get the collocation points from the C code, clean it up,
        find the function values at these points, apply the splitting function
        and make the C code do the spectral transfomation.

        Parameters
        ----------
        geometry : postcactus.grid_data.RegGeom
            The geometry at which the data should be read and interpolated
        quantiy : list
            A list of quantities which will be processed (redundant as of now)
        iterations : list
            A list of iteration the user want to process
        test : bool, optional
            Whether to make a test plot or not (default is False)
        split : bool optional
            Whether or not to apply the splitting function (default is True)
        scaling_factor : float, optional
            The scaling factor used to determine the distance at which
            the splitting factor will be used (default is 4.0)



        Notes
        -----

        - The code will now go through all quantities, so the quantity parameter is redundant and will be removed.
        - This code is as of now only applicable to equal-massed BHs. This will be change, and then there will be two scaling factors.
        """

        if self.nb_bodies > 1:
            possible_iterations, positions, radii = self.read_bbh_diag()

        if quantities == []:
            quantities = ["alp", "betax", "betay", "betaz",
                    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
                    "kxx", "kxy", "kxz", "kyy", "kyz", "kzz"]



        split = False if self.nb_bodies == 1 else split
        start_time = time.time()
        for it in iterations:
            print "[~] Starting with Iteration %s \n" %it


            if self.nb_bodies > 1:
                it_index = np.where(possible_iterations == it)

                if len(it_index) == 0:
                    continue

                pos1 = (positions[0,:,it_index])[0,0,:]
                pos2 = (positions[1,:,it_index])[0,0,:]
                radius1 = (radii[0,it_index])[0,0]
                radius2 = (radii[1,it_index])[0,0]

            else:
                pos1 = [0,0,0]
                pos2 = [0,0,0]

                radius1 = 0
                radius2 = 0

            
            if test:
                """
                if self.nb_bodies == 1:
                    self.make_test_plot("alp")
                elif self.nb_bodies == 2:
                    self.make_test_bbh_plot("gxx", pos2, pos1, radius2, radius1)
                """
                ETquantities.test_plot()



            for quantity in quantities:
                print "[~] Starting with Quantity %s" %quantity
                ETquantities.read(quantity)
                q = ETquantities#self.read_ET_quantity(quantity, g, it, dimentions=3, order=4)

                print "[+] Quantity Successfully Read from ET File"
                #BH1
                print "[~] Starting with Black Hole 1"


                values, flatten_values = self.get_values_at_coll_points(q,c_path=c_path,smooth=split, bh_pos=pos1, bh_rad=radius1,bh_pos2=pos2, bh_rad2=radius2, scaling_factor=scaling_factor)
                filename = "%s/%s_%s_body1.txt" %(result_path,quantity, it)
                self.write_flatten_values_to_file(flatten_values, it, 1, filename)


                #BH2
                if self.nb_bodies > 1:
                    print "[~] Now Black Hole 2"

                    values, flatten_values = self.get_values_at_coll_points(q,c_path=c_path,smooth=split, bh_pos=pos2, bh_rad=radius2,bh_pos2=pos1, bh_rad2=radius1, scaling_factor=scaling_factor)
                    filename = "%s/%s_%s_body2.txt" %(result_path,quantity, it)
                    self.write_flatten_values_to_file(flatten_values, it, 2, filename)

                print "\n INFO: Time used for %s: %.3f min.\n\n" %(quantities,(time.time()- start_time)/60.)

            if do_gyoto_converstion:
                print "[~] LORENE is Writing BH1 to GYOTO File"
                self.LORENE_read(filename, c_path=c_path, body=1, origin=pos1, it=it)

                if self.nb_bodies > 1:
                    print "[~] LORENE is Writing BH1 to GYOTO File"
                    self.LORENE_read(filename, c_path=c_path, body=2, origin=pos2, it=it)

            print "[+] Done with Iteration %s in %.3f min. \n" %(it, (time.time()- start_time)/60.)


        print "[+] Done in %.3f min!" %((time.time()- start_time)/60.)


    def _make_analytical_metric(self, analytical_function, c_path="../../C/", result_path="./", do_gyoto_converstion=True):
        quantities = ["alp", "betax", "betay", "betaz",
                    "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
                    "kxx", "kxy", "kxz", "kyy", "kyz", "kzz"]
        start_time = time.time()
        for quantity in quantities:
            
            print "[+] Function with metric component %s made" %quantity
            
            q = analytical_function(quantity)
            
            pos = np.array([0,0,0])


            values, flatten_values = self.get_values_at_coll_points(q,c_path=c_path,smooth=False, bh_pos=pos, bh_rad=1,bh_pos2=pos, bh_rad2=1, scaling_factor=0)
            filename = "%s/%s_%s_body1.txt" %(result_path,quantity, 0)
            self.write_flatten_values_to_file(flatten_values, 0, 1, filename)


            print "\n INFO: Time used: %.3f min.\n\n" %((time.time()- start_time)/60.)

        if do_gyoto_converstion:
            print "[~] LORENE is Writing BH1 to GYOTO File"
            self.LORENE_read(filename, c_path=c_path, body=1, origin=pos, it=0)

            


        print "[+] Done in %.3f min!" %((time.time()- start_time)/60.)


    def get_minkowski_component(self,quantity):
        if quantity in ["alp", "gxx", "gyy", "gzz"]:
            return lambda x : 1
        else:
            return lambda x : 0

    def _non_zero_schwarzschild_alp(self, x):
        r = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
        if abs(r) < 0.5*1e-6:
            r = 0.5*1e-6
        return (1-1/(2*r))/(1+1/(2*r)) #+ np.random.rand()/10**(4)
        
    def _non_zero_schwarzschild_gamma(self, x):
        r = sqrt(x[0]**2 + x[1]**2 + x[2]**2)
        if abs(r) < 0.5*1e-6:
            r = 0.5*1e-6
        return (1+1/(2*r))**4 #+ np.random.rand()/10**(4)


    def get_schwarzschild_isotropic_compinents(self, quantity):
        if quantity == "alp":
            return self._non_zero_schwarzschild_alp
        elif quantity in ["gxx", "gyy", "gzz"]:
            return self._non_zero_schwarzschild_gamma
        else:
            return lambda x : 0

        
    def make_minkowski(self, c_path="../../C/", result_path="./", do_gyoto_converstion=True):
        analytical_function = self.get_minkowski_component
        self._make_analytical_metric(analytical_function, c_path, result_path, do_gyoto_converstion)

    def make_schwarzschild_isotropic(self, c_path="../../C/", result_path="./", do_gyoto_converstion=True):
        analytical_function = self.get_schwarzschild_isotropic_compinents
        self._make_analytical_metric(analytical_function, c_path, result_path, do_gyoto_converstion)
    

    




if __name__=="__main__":
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/GW150914_28"
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/tov_ET_11"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/bh_3d"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_3D"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_large"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_large"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/schwarzschild_large"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_higherres"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires_center"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires_one"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires_center_hiphi"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_analytic"


    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr"
    folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/analytical_schwarz_cleaned_dx4"




    pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles"

    ### For making minkowski space
    inter = ETInterpolater("", 1)
    inter.xlim = [-10000,10000]
    inter.ylim = [-10000,10000]
    inter.zlim = [-10000,10000]
    #inter.make_minkowski(do_gyoto_converstion=False)
    inter.make_schwarzschild_isotropic(do_gyoto_converstion=False)
    exit()
    """
    """

    #quantity = "betay"
    #filename = "%s.txt" %quantity
    nb_bodies = 1

    inter = ETInterpolater(folder, nb_bodies)
    it = 0
    #g = inter.make_positive_geometry([-30,-30, -30], 30)
    """
    g0 = inter.make_positive_geometry([-10,-10, -10], 100)
    g15 = inter.make_positive_geometry([-20,-20, -20], 100)
    g1 = inter.make_positive_geometry([-100,-100, -100], 100)
    g2 = inter.make_positive_geometry([-300,-300, -300], 100)
    limits = [5,10, 50,250]
    et_q = ReadQuantities([g0,g15,g1,g2], it, folder, pickle_folder=pickle_folder, pickle=False, linear=True, limits=limits)
    """
    """
    """
    #g = inter.make_positive_geometry([-300,-300, -300], 300)
    g = None
    limits = [250]
    et_q = ReadQuantities([g], it, folder, pickle_folder=pickle_folder, pickle=False, linear=True, limits=limits)


    #et_q = ETQuantities(g, it, folder, pickle_folder=pickle_folder, pickle=False)
    #et_q = ETQuantities_gridInterpolator(g, it, folder, pickle_folder=pickle_folder, pickle=False)

    et_q.test_plot("gxx", binary=True)
    #exit()
    #inter.analyse_bbh(g, et_q, [it], quantities=["alp"], test=False)
    exit()

  
    exit()
    """
    index = d*(nr*nt*np) + k*nt*np + j*np + r
    print index, len(flatten_values), nz*nr*nt*np, (nz-1)*(nr*nt*np) + (nr-1)*nt*np + (nt-1)*np + np-1
    print(values[d][k][j][r])
    print(flatten_values[index])

    #inter.make_test_plot("alpha")

    #print(inter.get_coll_points()[0])
    """





"""
dx | RAM | Worked | Mean Error | Max Error


1  | 388.207 GByte | No |   |

1.5| 131.777 GByte | No |

1.875| 77.955 GByte |  Yes | 6.971825901943971e-05 | 0.0008034376108430052

2  | 65.708 GByte |  Yes | 7.450341308336404e-05 | 0.0009293562186072357

2.5| 37.045 GByte |  Yes | 0.00011172354624641227 | 0.0012555730776417917

3  | 23.927 GByte |  Yes | 0.00015018127776903771 | 0.0020727476781523535

4  | 12.728 GByte |  Yes | 0.0002152858966590524 | 0.003665644172405891


Large:

dx | RAM | Worked

2  | 155.498 GByte |  No

2.5| 83.485 GByte |  Yes

3  | 23.927 GByte |  Yes

4  | 12.728 GByte |  Yes
"""

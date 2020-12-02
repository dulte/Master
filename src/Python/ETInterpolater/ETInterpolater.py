import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import time

import os
import pickle

from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)


from postcactus.simdir import SimDir
from postcactus import visualize as viz
from postcactus import grid_data as gd

from scipy import ndimage


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

        corner = self.geo.x1()
        self.folder = pickle_folder
        self.filename = "et_quantities_%d_%d_%d" %(corner[0], corner[1], corner[2])

    def read(self, name):
        start_time = time.time()
        full_path = "%s/%s_%s" %(self.folder, name, self.filename)

        if self._check_pickle(name):
            print "[~] Pickle Found. Reading Pickle."
            self.load_pickle(name)

        else:
            print "[~] Pickle Not Found. Reading From ET Files."

            grid = self._read_quantity(name, self.geo, self.iteration)
            self.loaded_spline = grid.spline(order=3, mode="nearest")

        if self.pickle == True and not os.path.exists(full_path):
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

        n = 300
        end = 20
        q_inter = np.zeros((n,n))
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        for i in range(n):
            for j in range(n):
                input = np.array([abs(x[i]), abs(y[j]),0])

                q_inter[j,i] = self(input)


        plt.pcolormesh(x,y,q_inter)
        plt.colorbar()
        plt.show()

        plt.contour(x,y,q_inter)
        plt.show()


    def __call__(self, coords, output=None):


        q = self.loaded_spline
        if q is None:
            raise ValueError("[-] No Quantity Loaded")

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
        self.xlim = [0, abs(corner[0])]
        self.ylim = [0, abs(corner[1])]

        if len(corner) == 3:
            self.zlim = [0, abs(corner[2])]

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


    def get_values_at_coll_points(self, interpolated_quantity, smooth=True, bh_pos=[0,0,0], bh_rad=0,bh_pos2=[0,0,0], bh_rad2=0, scaling_factor=4.0, test=False):
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

        xx, yy, zz = self.get_coll_points(bh_pos)

        bbh_distance = np.linalg.norm(bh_pos-bh_pos2)


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
            if tp == "x":
                coord = self.xlim[0] - 5
            elif tp == "y":
                coord = self.ylim[0] - 5
            elif tp == "z":
                coord = self.zlim[0] - 5
            else:
                raise ValueError("%s is not an axis" %tp)
        elif coord == -np.inf:# or coord == -np.nan:
            #print coord, self.xlim[1]
            if tp == "x":
                coord = self.xlim[1] + 5
            elif tp == "y":
                coord = self.ylim[1] + 5
            elif tp == "z":
                coord = self.zlim[1] + 5
            else:
                raise ValueError("%s is not an axis" %tp)

        return self.desymmetrize_coord(coord)

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
                f.write("%.5f " %print_nb)

        print "[+] File Successfully Written"
        #np.savetxt(file, np.array(values))


    def LORENE_read(self, filename, c_path="../../C/" origin=[0,0,0], body=1, it=0):
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



    def get_coll_points(self,c_path="../../C/" origin=[0,0,0], body=1, it=0):
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
            raise IOError("Could not read LORENE C code!")
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





    def analyse_bbh(self, geometry, ETquantities, iterations, c_path="../../C/", result_path="", do_gyoto_converstion=True, quantities=[], test=False, split=True, scaling_factor=4.0):
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

                print "\n INFO: Time used: %.3f min.\n\n" %((time.time()- start_time)/60.)

            if do_gyoto_converstion:
                print "[~] LORENE is Writing BH1 to GYOTO File"
                self.LORENE_read(filename, c_path=c_path, body=1, origin=pos1, it=it)

                if self.nb_bodies > 1:
                    print "[~] LORENE is Writing BH1 to GYOTO File"
                    self.LORENE_read(filename, c_path=c_path, body=2, origin=pos2, it=it)

            print "[+] Done with Iteration %s in %.3f min. \n" %(it, (time.time()- start_time)/60.)


        print "[+] Done in %.3f min!" %((time.time()- start_time)/60.)







if __name__=="__main__":
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/GW150914_28"
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/tov_ET_11"
    folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_3D"

    pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles"

    
    
    

    
    quantity = "alp"
    filename = "%s.txt" %quantity
    inter = ETInterpolater(folder, 2)
    #g = inter.make_geometry([-50, -50, -50], 400)
    g = inter.make_positive_geometry([-20,-20, -20], 400)
    it = 0


    et_q = ETQuantities(g, it, folder, pickle_folder=pickle_folder, pickle=False)
    #et_q.read_all()
    #et_q.read("betax")
    et_q.test_plot("kxy")
    exit()
    inter.analyse_bbh(g, et_q, [it], quantities=["alp"], test=False)

    #q = inter.read_ET_quantity(quantity, g, it, dimentions=3, order=4)

    #inter.make_test_plot("gxx")
    exit()
    values, flatten_values = inter.get_values_at_coll_points(q, test=False)
    inter.write_flatten_values_to_file(flatten_values, it, 1, filename)
    inter.LORENE_read(filename)



    nz = 3
    nr = 8
    nt = 5
    np = 7



    d = 0
    k = 1
    j = 3
    r = 4
    print(values[d][k][j][r])
    """
    index = d*(nr*nt*np) + k*nt*np + j*np + r
    print index, len(flatten_values), nz*nr*nt*np, (nz-1)*(nr*nt*np) + (nr-1)*nt*np + (nt-1)*np + np-1
    print(values[d][k][j][r])
    print(flatten_values[index])

    #inter.make_test_plot("alpha")

    #print(inter.get_coll_points()[0])
    """

import numpy as np
import matplotlib.pyplot as plt
import h5py
import subprocess
import time

from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)


from postcactus.simdir import SimDir
from postcactus import visualize as viz
from postcactus import grid_data as gd


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
"""


class ETInterpolater:
    def __init__(self, dir, nb_bodies=2):
        self.sd = SimDir(dir)
        self.nb_bodies = nb_bodies

        self.dir = dir

        self.xlim = [0,0]
        self.ylim = [0,0]
        self.zlim = [0,0]


    def read_ET_quantity(self,quantity, geometry, iteration, dimentions=3, order=4):
        if dimentions == 2:
            grid = self.sd.grid.xy.read(quantity, iteration, geom=geometry, order=order)

        elif dimentions == 3:
            grid = self.sd.grid.xyz.read(quantity, iteration, geom=geometry, order=order)
        else:
            raise ValueError("Number of dimentions should be 2 or 3!")

        print "[+] %s Successfully Read" %quantity
        return grid


    def make_geometry(self, corner, n_pts):
        corner1 = [-i for i in corner]
        self.xlim = [corner1[0], corner[0]]
        self.ylim = [corner1[1], corner[1]]

        if len(corner) == 3:
            self.zlim = [corner1[2], corner[2]]
        geo = gd.RegGeom([n_pts]*len(corner), corner, x1=corner1)
        print "[+] Geomentry Successfully Made"
        return geo



    def make_test_plot(self, quantity):
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
        g = self.make_geometry([-50,-50], 400)
        q = self.read_ET_quantity(quantity, g, 0, dimentions=2)

        n = 600
        end = 10
        q_inter = np.zeros((n,n))
        x = np.linspace(-end, end, n)
        y = np.linspace(-end, end, n)

        th = np.linspace(0,2*np.pi,100)


        for i in range(n):
            for j in range(n):
                q_inter[j,i] = q(np.array([self.desymmetrize_coord(x[i]), self.desymmetrize_coord(y[j])]))
                q_inter[j,i] *= self.super_gaussian(x[i], y[j], 0, p1, r1)*self.super_gaussian(x[i], y[j], 0, p2, r2)


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


    def get_values_at_coll_points(self, interpolated_quantity,smooth=False, bh_pos=[0,0,0], bh_rad=0,bh_pos2=[0,0,0], bh_rad2=0, test=False):
        q = interpolated_quantity
        xx, yy, zz = self.get_coll_points(bh_pos)


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
                        x = self.bound_coord(xx[index][k][j][r], "x")
                        y = self.bound_coord(yy[index][k][j][r], "y")
                        z = self.bound_coord(zz[index][k][j][r], "z")



                        inter_q = q([x,y,z])
                        if smooth:
                            inter_q *= self.super_gaussian(x,y,z, bh_pos, bh_rad)
                            inter_q *= self.super_gaussian(x,y,z, bh_pos2, bh_rad2)

                        values[index][k][j][r] = inter_q

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

        if coord == np.inf or coord != coord:
            if tp == "x":
                coord = self.xlim[0] - 5
            elif tp == "y":
                coord = self.ylim[0] - 5
            elif tp == "z":
                coord = self.zlim[0] - 5
        elif coord == -np.inf:# or coord == -np.nan:
            #print coord, self.xlim[1]
            if tp == "x":
                coord = self.xlim[1] + 5
            elif tp == "y":
                coord = self.ylim[1] + 5
            elif tp == "z":
                coord = self.zlim[1] + 5

        return coord#self.desymmetrize_coord(coord)

    def desymmetrize_coord(self, coord):
        return abs(coord)



    def flatten_dict(self, d):
        flatten_values = []
        for index,domain in enumerate(d):
            for k in domain.keys():
                for j in domain[k].keys():
                    for r in range(len(domain[k][j])):
                        flatten_values.append(d[index][k][j][r])

        return flatten_values



    def write_values_to_file(self, values, file):
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
        values = [body, self.nb_bodies, it] + values
        with open(file, "w") as f:
            for i in values:
                print_nb = i
                if print_nb != print_nb:
                    print_nb = np.inf
                f.write("%.5f " %print_nb)

        print "[+] File Successfully Written"
        #np.savetxt(file, np.array(values))


    def LORENE_read(self, filename, origin=[0,0,0], body=1, it=0):
        p = subprocess.Popen("../C/get_points %s %s %s 1 %s %s" %(origin[0], origin[1], origin[2], body, it), stdout=subprocess.PIPE, shell=True)
        #p = subprocess.Popen("./get_points %s %s %s 1" %(origin[0], origin[1], origin[2]), stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        if p_status != 0:
            raise IOError("Could not read LORENE C code!")
        else:
            print "[+] LORENE Successfully Read and Saved the File"



    def get_coll_points(self, origin=[0,0,0], body=1, it=0):
        p = subprocess.Popen("../C/get_points %s %s %s 0 %s %s" %(origin[0], origin[1], origin[2], body, it), stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        if p_status != 0:
            raise IOError("Could not read LORENE C code!")
        else:
            print "[+] LORENE Successfully Gave the Coll Points"

        #print output
        x_dicts, y_dicts, z_dicts = self.clean_coll_points_xyz(output)
        return x_dicts, y_dicts, z_dicts



    def clean_coll_points_xyz(self,s):
        """
        Takes a string containing the output of the LORENE code and returns
        three dicts (x,y,z) giving the coord position of each coll point.

        Very bad practis is used in this function, so might need to be rewritten.
        """
        corrds = []
        index = 0
        for index, sub_s in enumerate(s):
                if sub_s == "+":
                    break
        s = s[index+1:]
        for index, sub_s in enumerate(s):
                if sub_s == "+":
                    break


        x_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        for index, sub_s in enumerate(s):
                if sub_s == "+":
                    break

        y_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        for index, sub_s in enumerate(s):
                if sub_s == "+":
                    break

        z_dicts = self.clean_coll_points(s[:index])
        s = s[index+1:]

        return x_dicts, y_dicts, z_dicts


    def clean_coll_points(self,s):
        """
        Takes the string containing the coll points for one coord (x,y or z),
        and returns a list with a dict for each domain.
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
        Takes a string for one domain of one coord and returns the dict for that domain.
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
        r = np.sqrt((bh_pos[0]-x)**2 + (bh_pos[1]-y)**2 + (bh_pos[2]-z)**2)/float(rad)
        return 1-I*np.exp(-2*r**n)



    def read_bbh_diag(self):
        itbh1, tbh1, xbh1, ybh1, zbh1, rbh1 = np.loadtxt("%s/output-0000/bbh3d/BH_diagnostics.ah1.gp" %self.dir,
                                                         usecols=(0,1,2,3,4,7),
                                                         unpack=True)

        itbh2, tbh2, xbh2, ybh2, zbh2, rbh2 = np.loadtxt("%s/output-0000/bbh3d/BH_diagnostics.ah2.gp" %self.dir,
                                                         usecols=(0,1,2,3,4,7),
                                                         unpack=True)


        return np.array(itbh1), np.array([[xbh1, ybh1, zbh1], [xbh2, ybh2, zbh2]]), np.array([rbh1, rbh2])





    def analyse_bbh(self, geometry, quantity, iterations, test=False):
        possible_iterations, positions, radii = self.read_bbh_diag()

        quantites = ["alp", "betax", "betay", "betaz",
                "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
                "kxx", "kxy", "kxz", "kyy", "kyz", "kzz"]
        start_time = time.time()
        for it in iterations:
            print "[~] Starting with Iteration %s \n" %it

            it_index = np.where(possible_iterations == it)

            if len(it_index) == 0:
                continue

            pos1 = (positions[0,:,it_index])[0,0,:]
            pos2 = (positions[1,:,it_index])[0,0,:]
            radius1 = (radii[0,it_index])[0,0]
            radius2 = (radii[1,it_index])[0,0]

            if test:
                self.make_test_bbh_plot(quantity, pos1, pos2, radius1, radius2)
            for quantity in quantites:
                print "[~] Starting with Quantity %s" %quantity
                q = self.read_ET_quantity(quantity, g, it, dimentions=3, order=4)

                print "[+] Quantity Successfully Read from ET File"
                #BH1
                print "[~] Starting with Black Hole 1"

                values, flatten_values = self.get_values_at_coll_points(q,smooth=True, bh_pos=pos1, bh_rad=radius1,bh_pos2=pos2, bh_rad2=radius2)
                filename = "../Python/%s_%s_body1.txt" %(quantity, it)
                self.write_flatten_values_to_file(flatten_values, it, 1, filename)


                #BH2
                print "[~] Now Black Hole 2"

                values, flatten_values = self.get_values_at_coll_points(q,smooth=True, bh_pos=pos2, bh_rad=radius2,bh_pos2=pos1, bh_rad2=radius1)
                filename = "../Python/%s_%s_body2.txt" %(quantity, it)
                self.write_flatten_values_to_file(flatten_values, it, 2, filename)

                print "\n INFO: Time used: %.3f min.\n\n" %((time.time()- start_time)/60.)


            print "[~] LORENE is Writing BH1 to GYOTO File"
            self.LORENE_read(filename, body=1, origin=pos1, it=it)

            print "[~] LORENE is Writing BH1 to GYOTO File"
            self.LORENE_read(filename, body=2, origin=pos2, it=it)

            print "[+] Done with Iteration %s in %.3f min. \n" %(it, (time.time()- start_time)/60.)


        print "[+] Done in %.3f min!" %((time.time()- start_time)/60.)







if __name__=="__main__":
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/GW150914_28"
    #folder = "/media/dulte/Seagate Expansion Drive/Storage/Master/ET_data/tov_ET_11"
    folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_3D"

    quantity = "kxx"
    filename = "%s.txt" %quantity
    inter = ETInterpolater(folder)
    g = inter.make_geometry([-40, -40, -40], 200)
    it = 0

    inter.analyse_bbh(g, quantity, [it], test=False)

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

#from ETInterpolater import ETInterpolater, ETQuantities, ETQuantities_gridInterpolator, ReadQuantities
from ET2G.ETInterpolater import ETInterpolater
from ET2G.ReadQuantities import ReadQuantities
from mpi4py import MPI
from time import sleep

# Sets up MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#folder = "/mn/stornext/d13/euclid/daniehei/simulations/bh_3d"
#pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles/smallTest"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_3D"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_large"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_large"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/schwarzschild_large"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_higherres"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires_center_2"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_analytic"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr_hires_one"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/analytical_schwarz_cleaned"



#folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
folder = "/mn/stornext/d13/euclid/daniehei/simulations/analytical_schwarz_cleaned_dx2"
pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles"

#quantities = ["gxx", "gyy"]

quantities = ["alp", "betax", "betay", "betaz",
        "gxx", "gxy", "gxz", "gyy", "gyz", "gzz",
        "kxx", "kxy", "kxz", "kyy", "kyz", "kzz"]


it = 0

nb_q_per_proc = len(quantities)//size
nb_rest = len(quantities)%size
nb_q_extra = 1 if rank < nb_rest else 0
nb_prev_extra = 1 if rank-1 < nb_rest else 0


rank_start = rank*(nb_q_per_proc)
rank_end = rank_start + nb_q_per_proc

rank_quantities = quantities[rank_start:rank_end]
if rank < nb_rest:
    rank_quantities.append(quantities[(size)*nb_q_per_proc+rank])

sleep(1)
sleep(rank)

print rank, rank_quantities

sleep(3)
nb_bodies = 1
linear = True
inter = ETInterpolater(folder, nb_bodies)
g0 = inter.make_positive_geometry([-10,-10, -10], 100)
g1 = inter.make_positive_geometry([-50,-50, -50], 100)
g2 = inter.make_positive_geometry([-200,-200, -200], 200)

it = 0

#limits = [5, 50, 400]
limits = [300]
g = None#inter.make_positive_geometry([-500,-500, -500], 250)

"""
For best image of kerr use [-200,-200, -200], 100
"""


#et_q = ReadQuantities([g0,g1,g2], it, folder, pickle_folder=pickle_folder, pickle=False, linear=linear,limits=limits)
et_q = ReadQuantities([g], it, folder, pickle_folder=pickle_folder, pickle=False, linear=linear,limits=limits)
#et_q = ETQuantities(g, it, folder, pickle_folder=pickle_folder, pickle=False)
#et_q = ETQuantities_gridInterpolator(g, it, folder, pickle_folder=pickle_folder, pickle=False)


smooth = True

inter.analyse_bbh(g, et_q, [it],quantities=rank_quantities, test=False, do_gyoto_converstion=False, split=smooth)

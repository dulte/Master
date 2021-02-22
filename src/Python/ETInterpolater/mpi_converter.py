from ETInterpolater import ETInterpolater, ETQuantities, ETQuantities_gridInterpolator
from mpi4py import MPI
from time import sleep

# Sets up MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
#pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles/smallTest"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/tov_3D"
#folder = "/mn/stornext/d13/euclid/daniehei/simulations/kerr"
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
nb_bodies = 2
inter = ETInterpolater(folder, nb_bodies)
g = inter.make_positive_geometry([-80,-80, -80], 400)
#g = inter.make_positive_geometry([-200,-200, -200], 800)

"""
For best image of kerr use [-200,-200, -200], 100
"""

et_q = ETQuantities(g, it, folder, pickle_folder=pickle_folder, pickle=False)
#et_q = ETQuantities_gridInterpolator(g, it, folder, pickle_folder=pickle_folder, pickle=False)


smooth = True

inter.analyse_bbh(g, et_q, [it],quantities=rank_quantities, test=False, do_gyoto_converstion=False, split=smooth)

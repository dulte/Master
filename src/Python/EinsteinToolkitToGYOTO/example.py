from ET2G.ETInterpolater import ETInterpolater
from ET2G.ReadQuantities import ReadQuantities

if __name__=="__main__":
    #######################################################
    ####      Example of making analytical metrics     ####                                       
    #######################################################
    
    inter = ETInterpolater("", 1)
    inter.xlim = [-10000,10000]
    inter.ylim = [-10000,10000]
    inter.zlim = [-10000,10000]
    #inter.make_minkowski(do_gyoto_converstion=False)
    inter.make_schwarzschild_isotropic(do_gyoto_converstion=False)
    exit()
    """
    """

    #######################################################
    ##  Example of reading and converting simulated data ##                                       
    #######################################################

    #folder = "/mn/stornext/d13/euclid/daniehei/simulations/bbh_3D"
    folder = "/mn/stornext/d13/euclid/daniehei/simulations/analytical_schwarz_cleaned_dx2"
    #pickle_folder = "/mn/stornext/d13/euclid/daniehei/ETConverter/spline_pickles"



    nb_bodies = 1
    inter = ETInterpolater(folder, nb_bodies)
    it = 0
    
    g = None
    
    et_q = ReadQuantities([g], it, folder, pickle_folder=pickle_folder, pickle=False, linear=True, limits=limits)


    et_q.test_plot("gxx", binary=True)
    exit()
    
    
    
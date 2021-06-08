import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from math import sqrt
import seaborn as sns

sns.set_style("darkgrid")
plt.rcParams.update({"font.size": 30})


def get_norm(filename):

    norm_gyoto = []

    with open(filename) as f:
        prev_r = 10000
        photon_nr = 0
        closet_r = 100000
        for line in f:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0] == "Norm:":
                
                n = float(words[1])
                r = float(words[2])
                p = float(words[3])
                t = float(words[4])
                tdot = float(words[5])


                if r > 200 and prev_r < 100:
                    
                    closet_r = r
                    photon_nr += 1
                    prev_r = r

                if r < prev_r:
                    closet_r = r
            
                prev_r = r
                s = [n, r, p, t, tdot, photon_nr, closet_r]
                
                norm_gyoto.append(s)

    return np.array(norm_gyoto)


def get_mean(norm, scale=10,use_tdot=False):

    result = {i/float(scale):[] for i in range(1,250*scale + 1)}



    for i in range(len(norm)):
        r = norm[i][1]
        n = norm[i][0]
        tdot = norm[i][4]
        r_key = round(r*scale)/float(scale)
        if use_tdot:
            n /= tdot**2
        result[r_key].append(n)


    mean_result = [np.mean(np.abs(result[key])) if len(result[key])>0 else 0 for key in result.keys() ]

    return result.keys(), np.array(mean_result)



#filename_gyoto = "norm_case5-nostar.txt"
"""
case = 6

norm_gyoto = get_norm("norm_case%s-nostar.txt" %case)
norm_analytic = get_norm("norm_case%s-star-20.txt" %case)
norm_lorene = get_norm("norm_case%s-star-5.txt" %case)
"""

norm_gyoto = get_norm("norm_case3-nostar.txt")
norm_analytic = get_norm("norm_case7-nostar.txt")
norm_lorene = get_norm("norm_case8-nostar.txt")

rs_ET = np.load("r_ET_data.npy")
mean_diff = np.load("mean_diff_ET_data.npy")

rs, mean_result = get_mean(norm_gyoto,use_tdot=True)
rs_a, mean_result_a = get_mean(norm_analytic,use_tdot=True)
rs_l, mean_result_l = get_mean(norm_lorene,use_tdot=True)



photon_nr = [30*15 + i for i in range(7)]
single_photon = norm_gyoto[np.where(np.logical_and(norm_gyoto[:,5] >= photon_nr[0], norm_gyoto[:,5] <= photon_nr[-1]))]
single_photon_a = norm_analytic[np.where(np.logical_and(norm_analytic[:,5] >= photon_nr[0], norm_analytic[:,5] <= photon_nr[-1]))]

fig = plt.figure()
ax = fig.add_subplot(projection='polar')
polar_cmap = plt.cm.get_cmap("rocket").reversed()

#sc = ax.scatter(norm_gyoto[:,3], norm_gyoto[:,1], c=np.log(np.abs(norm_gyoto[:,0])),marker=".",cmap='hsv')
sc = ax.scatter(single_photon_a[:,3], single_photon_a[:,1], c=np.abs(single_photon_a[:,0]),marker="+",norm=LogNorm(vmin=10e-8,vmax=10e-1),cmap=polar_cmap,label="LORENE")
sc = ax.scatter(single_photon[:,3], single_photon[:,1], c=np.abs(single_photon[:,0]),marker=".",norm=LogNorm(vmin=10e-8,vmax=10e-1),cmap=polar_cmap,label="Simulated")
ax.set_rmax(30)
circle = plt.Circle((0.0, 0.0), 0.5, transform=ax.transData._b, color="black")
ax.add_artist(circle)
plt.legend()
plt.colorbar(sc)
plt.show()


normalized_gyoto = abs(norm_gyoto[:,0])/norm_gyoto[:,4]**2
normalized_analytic = abs(norm_analytic[:,0])/norm_analytic[:,4]**2
normalized_lorene = abs(norm_lorene[:,0])/norm_lorene[:,4]**2

plt.semilogy(norm_gyoto[:,1], normalized_gyoto,"+",label="Simulated Case 1",alpha=0.5)
plt.semilogy(norm_analytic[:,1], normalized_analytic,"x",label="Simulated Case 2",alpha=0.3)
plt.semilogy(norm_lorene[:,1], normalized_lorene,"+",label="LORENE",alpha=0.3)
"""
plt.semilogy(norm_gyoto[:,1], normalized_gyoto,"+",label="No Star",alpha=0.5)
plt.semilogy(norm_analytic[:,1], normalized_analytic,"x",label="Star RMax=20",alpha=0.3)
plt.semilogy(norm_lorene[:,1], normalized_lorene,"+",label="Star RMax=5",alpha=0.3)
"""

plt.legend()
plt.title(r"Photon Drift/$\dot{t}^2$")
plt.ylabel(r"Photon Drift/$\dot{t}^2$")
plt.xlabel("r")
plt.show()



infalling_photons = np.where(norm_gyoto[:,1]<=norm_gyoto[:,6])
outgoing_photons = np.where(norm_gyoto[:,1]>norm_gyoto[:,6])


fig, ax = plt.subplots(figsize=[5, 4])
axins = ax.inset_axes([0.15, 0.05, 0.4, 0.4])
ax.semilogy(norm_gyoto[outgoing_photons,1][0], normalized_gyoto[outgoing_photons],".",label="Outgoing Photon")
ax.semilogy(norm_gyoto[infalling_photons,1][0], normalized_gyoto[infalling_photons],".",label="Infalling Photons",alpha=0.5)
axins.semilogy(norm_gyoto[outgoing_photons,1][0], normalized_gyoto[outgoing_photons],".",label="Outgoing Photon")
axins.semilogy(norm_gyoto[infalling_photons,1][0], normalized_gyoto[infalling_photons],".",label="Infalling Photons",alpha=0.5)
x1, x2, y1, y2 = 0, 30, 1e-8, 1e-4
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels('')
axins.set_yticklabels('')
plt.legend()
plt.title(r"Difference in Norm of Ingoing and Outfalling Photons")
plt.ylabel(r"Photon Drift/$\dot{t}^2$")
plt.xlabel("r")
ax.indicate_inset_zoom(axins)
plt.show()





plt.scatter(norm_gyoto[:,1], normalized_gyoto, c=norm_gyoto[:,6],marker="+",cmap="rocket",alpha=0.5,norm=Normalize(vmin=0,vmax=20))
plt.yscale("log")
plt.colorbar(label="Nearest r")
plt.xlim([0,30])
plt.title(r"Photon Drift/$\dot{t}^2$")
plt.ylabel(r"Photon Drift/$\dot{t}^2$")
plt.xlabel("r")
plt.title("The Norm as a Function of r and the Nearest r")
plt.show()



plt.semilogy(rs, mean_result, ".",label="Simulated Case 1")
plt.semilogy(rs_a, mean_result_a, ".",label="Simulated Case 2",alpha=0.5)
plt.semilogy(rs_l, mean_result_l, ".",label="LORENE",alpha=0.2)
"""
plt.semilogy(rs, mean_result, ".",label="No Star")
plt.semilogy(rs_a, mean_result_a, ".",label="Star RMax=20",alpha=0.5)
plt.semilogy(rs_l, mean_result_l, ".",label="Star RMax=5",alpha=0.2)
"""

plt.legend()
plt.title("Photon Drift/$\dot{t}^2$ over Bins size 0.1")
plt.ylabel("Photon Drift/$\dot{t}^2$")
plt.xlabel("r")
#plt.semilogy(rs_ET, mean_diff*(max(mean_result)/np.max(mean_diff)))
plt.show()

plt.semilogy(norm_gyoto[:,1], norm_gyoto[:,4],"+",label="Simulated",alpha=0.5)
plt.semilogy(norm_analytic[:,1], norm_analytic[:,4],"x",label="Analytic",alpha=0.3)
plt.legend()
plt.title(r"Time derivative of $t$: $\dot{t}^2$")
plt.ylabel(r"$\dot{t}^2$")
plt.xlabel("r")
plt.show()
"""
"""

plt.semilogy(rs, np.abs(mean_result/mean_result_l), ".")
plt.show()
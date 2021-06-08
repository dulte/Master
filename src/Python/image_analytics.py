import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from astropy.io import fits
import seaborn as sns

sns.set_style("dark")




def get_image(name):
    hdul = fits.open(name)
    return np.array(hdul[0].data)[0]


def plot_individual(image_a, image_s_a, metric_type, metric_type_2, compare_lorene=False):

    #Errors
    print(np.sum(np.abs(image_s_a-image_a))/np.sum(image_a))
    print(np.sum(np.abs(image_s_a-image_a))/900)
    print(np.sum(np.abs(image_s_a-image_a)))

    #Individual Plots:
    plt.rcParams.update({"font.size": 30})

    if not  compare_lorene:
        plt.imshow(image_a)
        plt.title(r"Intesity  $\mathcal{I}(x,y)$ for " + metric_type)
        plt.xlabel("Pixel x")
        plt.ylabel("Pixel y")
        plt.colorbar(extend="max",label=r"Intensity")
        plt.show()


        plt.imshow(image_s_a)
        plt.title(r"Intesity  $\mathcal{I}(x,y)$ for " + metric_type_2)
        plt.xlabel("Pixel x")
        plt.ylabel("Pixel y")
        plt.colorbar(extend="max",label=r"Intensity")
        plt.show()

    plt.imshow(np.abs(image_s_a-image_a), norm=LogNorm())
    #plt.title(r"Difference in Intesity  $\mathcal{I}_{Simulated}-\mathcal{I}_{Analytic}$")
    if compare_lorene:
        plt.title(r"Difference in Intesity  $\mathcal{I}-\mathcal{I}_{LORENE}$")
    else:
        plt.title(r"Difference in Intesity  RMax=20 and RMax=5")
    plt.xlabel("Pixel x")
    plt.ylabel("Pixel y")
    plt.colorbar(extend="max",label=r"Difference in Intensity")
    plt.show()

case = 5



print("Comparing Different RMmax")


image_a = get_image("case%s-star-5.fits" %case)
image_s_a = get_image("case%s-star-20.fits" %case)

metric_type = "RMax=5"
metric_type_2 = "RMax=20"

plot_individual(image_a, image_s_a, metric_type, metric_type_2, compare_lorene=False)

print("Comparing with LORENE")

image_a = get_image("output_page-thorn_lorene.fits" %case)
image_s_a = get_image("output_page-thorn_lorene_large.fits")

metric_type = "Test Metic"
metric_type_2 = "LORENE Metic"

plot_individual(image_a, image_s_a, metric_type, metric_type_2, compare_lorene=True)





























"""
#All in one plot
gs = gridspec.GridSpec(1, 3)

axs = [
    plt.subplot(gs[0,0]),
    plt.subplot(gs[0,1]),
    plt.subplot(gs[0,2]),
]



im1 = axs[0].imshow(image_a)
axs[0].set_title(r"Intesity  $\mathcal{I}(x,y)$ for " + metric_type)
axs[0].set(xlabel="Pixel x", ylabel="Pixel x")
#axs[0,0].xlabel("Pixel x")
#axs[0,0].ylabel("Pixel y")
plt.colorbar(im1, ax=axs[0], extend="max",label=r"Intensity")
#plt.show()


im2 = axs[1].imshow(image_s_a)
axs[1].set_title(r"Intesity  $\mathcal{I}(x,y)$ for " + metric_type_2)
axs[1].set(xlabel="Pixel x", ylabel="Pixel y")
#axs[0,1].xlabel("Pixel x")
#axs[0,1].ylabel("Pixel y")
plt.colorbar(im1, ax=axs[1], extend="max",label=r"Intensity")
#plt.show()

im3 = axs[2].imshow(np.abs(image_s_a-image_a), norm=LogNorm())
axs[2].set_title(r"Difference in Intesity  $\mathcal{I}_{Simulated}-\mathcal{I}_{Analytic}$")
axs[2].set(xlabel="Pixel x", ylabel="Pixel y")
#axs[1,0].xlabel("Pixel x")
#axs[1,0].ylabel("Pixel y")
plt.colorbar(im3, ax=axs[2], extend="max",label=r"Difference in Intensity")
plt.show()
"""


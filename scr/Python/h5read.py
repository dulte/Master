import h5py

f = h5py.File("trK.0.xy.h5", "r")
for i in list(f.keys()):
    print(i)


dset = f["ML_BSSN::trK it=0 tl=0 m=0 rl=0"]
print(dset[100,1])

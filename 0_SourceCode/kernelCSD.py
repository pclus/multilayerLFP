
import numpy as np
import matplotlib.pyplot as plt
from kcsd import KCSD2D
from kcsd import oKCSD2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---- Functions ----

# --
# Read the bindary file "pre.bin"
def read_data_chunk(t0,tf):
    dt=1/2500.0
    m0=int(np.floor(t0/dt)-1)
    mf=int(np.ceil(tf/dt))
    mm = mf-m0
    n=384
    m = 2250000; 
    data=np.zeros((n,mm))
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/filtered_pre.bin", "rb") as fin:
        data[0,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='',offset=m0)
        for i in range(1,n):
            data[i,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='', offset=8*(m-mm))
    return data;
# --

# --
# modifyied from the KCSD tutorial.
def do_kcsd(ele_pos, pots):
    h = 10.0     # thickness
    sigma = 1.0 # S/m
    length = pots.shape[1] # size of time dimensions
    pots = pots.reshape((len(ele_pos), length)) 
    k = KCSD2D(ele_pos, pots, h=h, sigma=sigma,                                                                                                                                                       
               xmin=-200.0, xmax=200.0,
               ymin=-60.0, ymax=3900.0,
               n_src_init=2000, src_type='gauss', R_init=40,lambd=1e-7, gdx=5.0,gdy=5.0) # gdx and gdy are the spatial resolutions

    return k
# --

# ---
# compute error:
def kcsd_error(k,ele_x,ele_y):
    ownk = oKCSD2D(ele_pos, pots, h=k.h, sigma=k.sigma,                                                                                                                                                       
                xmin=k.xmin, xmax=k.xmax,
                ymin=k.ymin, ymax=k.ymax,
                n_src_init=k.n_src_init, src_type=k.src_type, R_init=k.R_init,lambd=k.lambd,own_est=np.array((ele_x,ele_y)),own_src=(k.src_x,k.src_y)) 

    est_pots = ownk.values('POT')
    error = np.linalg.norm(k.pots-est_pots)**2 + k.lambd*np.linalg.norm(est_pots);
    return ownk,error;
# ---

# # --
# # modifyied from the KCSD tutorial
# # can be used with: ax=make_plot(k.estm_x, k.estm_y, est_csd[:, :, 0], cmap=cm.bwr) 
# def make_plot(xx, yy, zz, title='CSD', cmap=cm.bwr):
#     fig = plt.figure(figsize=(7, 7))
#     ax = plt.subplot(111)
#     t_max = np.max(np.abs(zz))
#     levels = np.linspace(-1 * t_max, t_max, 32)
#     im = ax.contourf(xx, yy, zz, levels=levels, cmap=cmap)
#     ax.set_xlabel('X (mm)')
#     ax.set_ylabel('Y (mm)')
#     ax.set_title(title)
#     ticks = np.linspace(-1 * t_max, t_max, 3, endpoint=True)
#     plt.colorbar(im, orientation='horizontal', format='%.2f', ticks=ticks)
#     return ax
# # -- 

# ---- Workflow ----

# -- 
# Define electrode locations (centers of the 12x12 squares)
# x_pos = np.array([29, 61, 13, 45 ])   # 0 is the left boundary of the probe
x_pos = np.array([-6, 26, -22, 10 ])    # 0 is the center of the probe
ele_x = np.tile(x_pos,96)
ele_y = [i for i in range(6,3840, 20) for j in range(2)]
ele_pos = np.column_stack((ele_x, ele_y))
# --

# --
# Save and plot locations:
# np.savetxt('electrode_locations.dat',ele_pos)
# p1=plt.scatter(ele_pos[:, 0], ele_pos[:, 1], 10, c='k')
# plt.show()
# --


# -- 
# Read s seconds, compute, and plot the kCSD
s = 1.0
pots =read_data_chunk(0.0004,10.0)
k = do_kcsd(ele_pos, pots)
est_csd = k.values('CSD')
redk, err = kcsd_error(k,ele_x,ele_y)

# # --
# f1=plt.figure(1)
# plt.imshow(np.transpose(est_csd[:,::-1,0]),cmap=cm.bwr,aspect='auto') 
# plt.show()
# # --

# --
# save data to plot in gnuplot
# a=k.estm_x.reshape(k.estm_x.size,1)
# b=k.estm_y.reshape(k.estm_y.size,1)
# c=est_csd.reshape(est_csd.size,1)
# np.savetxt("temp.dat", np.array((a,b,c))[:,:,0].transpose())

# aux=ownk.values('CSD')
# a=ownk.estm_x.reshape(ownk.estm_x.size,1)
# b=ownk.estm_y.reshape(ownk.estm_y.size,1)
# c=aux.reshape(aux.size,1)
# np.savetxt("temp2.dat", np.array((a,b,c))[:,:,0].transpose())
# --


# --
# animation, does not work in VS, use terminal
# for i in range(2500):
#     plt.imshow(np.transpose(est_csd[:,::-1,i]),cmap=cm.bwr,aspect='auto') 
#     plt.pause(0.001)
#
# plt.show()
# --

# --
""""
Since we cannot load the whole time series we'll have to take chuncks of 100 miliseconds,
store them in different binaries, and then reorder the chunkcs in a unique big fatty binary file
"""" 



# est_csd = k.values('POT')
# np.savetxt('temp.dat',est_csd[:,:,0])


# # f2=plt.figure(2)
# # ax=make_plot(k.estm_x, k.estm_y, est_csd[:, :, 0], 
# #           title='Estimated CSD without CV', cmap=cm.bwr) # First time point
# # #f2.show()

# # check the preprint. it works, but it does not compute the error function, but a modified version
# k.cross_validate(lambdas=np.logspace(-9, -1, num=9), Rs=np.linspace(10, 100, num=10))

# CHECK SOURCE LOCATIONS:...
p1=plt.scatter(k.src_x, k.src_y, 10, c='k')
plt.show()


# est_csd = k.values('CSD')

# plt.imshow(p)
# plt.show()

# plt.plot(np.transpose(k.errs[:,:])); plt.show()

# np.savetxt('cv.dat',k.errs)
# np.savetxt('kCSD.dat',est_csd[:,:,0])



# # position and distance of estimated sources from the electrodes:

# plt.scatter(k.src_x,k.src_y);
# plt.imshow(k.src_estm_dists)
# plt.show()

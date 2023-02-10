
import numpy as np
import matplotlib.pyplot as plt

from kcsd import KCSD2D

# -- 
# Define electrode locations (centers of the 12x12 squares)
x_pos = np.array([29, 61, 13, 45 ])
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

# -- read a single data point [THIS SHOUT BE REMOVED SOON]
def read_data_point(j):
   n= 384
   m = 2250000; 
   data=np.zeros(n)
   with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/pre.bin", "rb") as fin:
       data[0]=np.fromfile(fin, dtype=np.double, count=1, sep='',offset=8*(j-1))
       for i in range(1,n):
           data[i]=np.fromfile(fin, dtype=np.double, count=1, sep='', offset=8*(m-1))
   return data;
# --

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
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/pre.bin", "rb") as fin:
        data[0,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='',offset=m0)
        for i in range(1,n):
            data[i,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='', offset=8*(m-mm))
    return data;
# --

pots =read_data_chunk(0.0004,100.0)

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def make_plot(xx, yy, zz, title='True CSD', cmap=cm.bwr):
    fig = plt.figure(figsize=(7, 7))
    ax = plt.subplot(111)
    # ax.set_aspect('equal')
    t_max = np.max(np.abs(zz))
    levels = np.linspace(-1 * t_max, t_max, 32)
    im = ax.contourf(xx, yy, zz, levels=levels, cmap=cmap)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_title(title)
    ticks = np.linspace(-1 * t_max, t_max, 3, endpoint=True)
    plt.colorbar(im, orientation='horizontal', format='%.2f', ticks=ticks)
    return ax


def do_kcsd(ele_pos, pots):
    h = 10.  # thickness
    sigma = 1.0 # S/m
    # pots = pots.reshape((len(ele_pos), 1)) # one time point 
    # k = KCSD2D(ele_pos, pots, h=h, sigma=sigma,                                                                                                                                                       
    #            xmin=-30.0, xmax=100.0,
    #            ymin=-60.0, ymax=3900.0,
    #            n_src_init=1000, src_type='gauss', R_init=40,lambd=1e-5, gdx=5.0,gdy=5.0) # gdx and gdy are the spatial resolutions

    length = pots.shape[1]
    pots = pots.reshape((len(ele_pos), length)) # 2500 time points
    k = KCSD2D(ele_pos, pots, h=h, sigma=sigma,                                                                                                                                                       
               xmin=-30.0, xmax=100.0,
               ymin=-60.0, ymax=3900.0,
               n_src_init=1000, src_type='gauss', R_init=40,lambd=1e-5, gdx=5.0,gdy=5.0) # gdx and gdy are the spatial resolutions
    return k


k = do_kcsd(ele_pos, pots)
est_csd = k.values('CSD')

# est_csd = k.values('POT')
# np.savetxt('temp.dat',est_csd[:,:,0])

# --
f1=plt.figure(1)
plt.imshow(np.transpose(est_csd[:,::-1,122]),cmap=cm.bwr,aspect='auto') 
plt.show()
# --

# --
# animation, does not work in VS, use terminal
# for i in range(2500):
#     plt.imshow(np.transpose(est_csd[:,::-1,i]),cmap=cm.bwr,aspect='auto') 
#     plt.pause(0.001)
#
# plt.show()
# --

ax=make_plot(k.estm_x, k.estm_y, est_csd[:, :, 0], 
          title='Estimated CSD without CV', cmap=cm.bwr) # First time point
plt.show()


# # f2=plt.figure(2)
# # ax=make_plot(k.estm_x, k.estm_y, est_csd[:, :, 0], 
# #           title='Estimated CSD without CV', cmap=cm.bwr) # First time point
# # #f2.show()

# # check the preprint. it works, but it does not compute the error function, but a modified version
# k.cross_validate(lambdas=np.logspace(-7, -1, num=7), Rs=np.linspace(10, 100, num=10))

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

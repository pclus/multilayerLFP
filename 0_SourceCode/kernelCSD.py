
import numpy as np
import matplotlib.pyplot as plt

from kcsd import KCSD2D


x_pos = np.array([29, 61, 13, 45 ])
ele_x = np.tile(x_pos,96)
ele_y = [i for i in range(6,3840, 20) for j in range(2)]
ele_pos = np.column_stack((ele_x, ele_y))
# p1=plt.scatter(ele_pos[:, 0], ele_pos[:, 1], 10, c='k')
# plt.show()
#

n=384
m = 2250000; 
data=np.zeros(n)
with open("Documents/Data/UPO-tACs/1_Raw/pre.bin", "rb") as fin:
    data[0]=np.fromfile(fin, dtype=np.double, count=1, sep='')
    for i in range(1,n):
        data[i]=np.fromfile(fin, dtype=np.double, count=1, sep='', offset=8*(m-1))


    
pots=data

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
    pots = pots.reshape((len(ele_pos), 1)) # first time point 
    k = KCSD2D(ele_pos, pots, h=h, sigma=sigma,                                                                                                                                                       
               xmin=0.0, xmax=70.0,
               ymin=0.0, ymax=3840.0,
               n_src_init=1000, src_type='gauss', R_init=25)
    return k

k = do_kcsd(ele_pos, pots)
est_csd = k.values('CSD')

ax=make_plot(k.estm_x, k.estm_y, est_csd[:, :, 0], 
          title='Estimated CSD without CV', cmap=cm.bwr) # First time point
# plt.show()

k.cross_validate(lambdas=None, Rs=np.linspace(2, 100.0, 48))

est_csd = k.values('CSD')
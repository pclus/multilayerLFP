
import numpy as np
import matplotlib.pyplot as plt
from kcsd import KCSD2D
from kcsd import oKCSD2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---- Functions ----
class kcsd_opts:
    def __init__(self, h=1.0,sigma=1.0,xmin=-200.0,xmax=200.0,ymin=-60.0,
                 ymax=3900.0,n_src_init=2000,src_type='gauss',R_init=40.0,lambd=1e-6,gdx=5.0,gdy=5.0):
        self.h=h
        self.sigma=sigma
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.n_src_init=n_src_init
        self.src_type=src_type
        self.R_init=R_init
        self.lambd=lambd
        self.gdx=gdx
        self.gdy=gdy

opts=kcsd_opts()
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
def do_kcsd(ele_pos, pots,opts):
    h = 10.0     # thickness
    sigma = 1.0 # S/m
    length = pots.shape[1] # size of time dimensions
    pots = pots.reshape((len(ele_pos), length)) 
    k = KCSD2D(ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
               xmin=opts.xmin, xmax=opts.xmax,
               ymin=opts.ymin, ymax=opts.ymax,
               n_src_init=opts.n_src_init, src_type=opts.src_type,
               R_init=opts.R_init,lambd=opts.lambd,
               gdx=opts.gdx,gdy=opts.gdy)

    return k
# --

# ---
# compute error:
def kcsd_error(ele_pos,pots,ele_x,ele_y,opts):
    ownk = oKCSD2D(ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
                xmin=opts.xmin, xmax=opts.xmax,
                ymin=opts.ymin, ymax=opts.ymax,
                n_src_init=opts.n_src_init, src_type=opts.src_type,
                R_init=opts.R_init,lambd=opts.lambd,
                own_est=np.array((ele_x,ele_y)),own_src=(opts.src_x,opts.src_y)) 

    est_pots = ownk.values('POT')
    error = np.linalg.norm(k.pots-est_pots)**2 + opts.lambd*np.linalg.norm(est_pots);
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
pots =read_data_chunk(0.0004,0.01)
k = do_kcsd(ele_pos, pots,opts)
est_csd = k.values('CSD')
opts.src_x=k.src_x;
opts.src_y=k.src_y;
redk, err = kcsd_error(ele_pos,pots,ele_x,ele_y,opts)
# --

# --
# CV to obtain reasonable values of R and lambda.
# It does not compute the error function, but a modified version, a bit heavy for long time series
k.cross_validate(lambdas=np.logspace(-9, -1, num=9), Rs=np.linspace(10, 100, num=10))
# Result: R,lambda : 60.0 1e-6
# --

# --
f1=plt.figure(1)
plt.imshow(np.transpose(est_csd[:,::-1,0]),cmap=cm.bwr,aspect='auto') 
plt.show()
# --

# --
def kcsd_validation(Rs,lambdas,k):
    kcsd_errors = np.zeros([Rs.size,lambdas.size])
    for i,R in enumerate(Rs):
        for j,lambd in enumerate(lambdas):
            kred, err = kcsd_error(k,ele_x,ele_y,R,lambd)
            kcsd_errors[i,j] = err;
    return kcsd_errors;
# --

lambdas=np.logspace(-9, -1, num=9);
Rs=np.linspace(10, 100, num=10);
errors=kcsd_validation(Rs,lambdas,ele_pos,pots)
plt.imshow(errors); plt.show();
# provides minimal values for R=40 and lambda=1e-9 ...but one should take into account the effect of lambda on the calculation of the error
# --

def export_at_electrodes(ele_pos,pots,ele_x,ele_y,ops):
    pots =read_data_chunk(0.0004,900.0)
    redk = oKCSD2D(ele_pos, pots, h=ops.h, sigma=ops.sigma,                                                                                                                                                       
            xmin=ops.xmin, xmax=ops.xmax,
            ymin=ops.ymin, ymax=ops.ymax,
            n_src_init=ops.n_src_init, src_type=ops.src_type, 
            R_init=ops.R_init,lambd=ops.lambd,own_est=np.array((ele_x,ele_y)),own_src=(ops.src_x,ops.src_y)) 
    est_csd = redk.values('CSD')
    # est_pot = redk.values('POT')
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/kCSD_electrodes_pre.bin", "wb") as fout:
        est_csd.tofile(fout,"");
    return est_csd;

def export_at_centers(ele_pos,pots,ops):
    pots =read_data_chunk(0.0004,900.0)
    loc_y=np.arange(opts.ymin,opts.ymax,ops.gdy)
    loc_x=np.zeros(loc_y.size)
    redk = oKCSD2D(ele_pos, pots, h=ops.h, sigma=ops.sigma,                                                                                                                                                       
            xmin=ops.xmin, xmax=ops.xmax,
            ymin=ops.ymin, ymax=ops.ymax,
            n_src_init=ops.n_src_init, src_type=ops.src_type, 
            R_init=ops.R_init,lambd=ops.lambd,own_est=np.array((loc_x,loc_y)),own_src=(ops.src_x,ops.src_y)) 
    est_csd = redk.values('CSD')
    # est_pot = redk.values('POT')
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/kCSD_centers_pre.bin", "wb") as fout:
        est_csd.tofile(fout,"");
    return est_csd;

csd_electro=export_at_electrodes(ele_pos,pots,ele_x,ele_y,opts)
csd_centers=export_at_centers(ele_pos,pots,opts)

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



# CHECK SOURCE LOCATIONS:...
# p1=plt.scatter(k.src_x, k.src_y, 10, c='k')
# plt.show()


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

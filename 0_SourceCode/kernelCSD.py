
import numpy as np
import matplotlib.pyplot as plt
from kcsd import KCSD2D
from kcsd import oKCSD2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class kcsd_opts:
    def __init__(self, h=10.0,sigma=1.0,xmin=-200.0,xmax=200.0,ymin=-60.0,
                 ymax=3900.0,n_src_init=2000,src_type='gauss',R_init=60.0,lambd=1e-6,gdx=5.0,gdy=5.0):
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
        # -- 
        # Define electrode locations (centers of the 12x12 squares)
        self.x_pos = np.array([-6, 26, -22, 10 ])    # x=0 is the center of the probe
        self.ele_x = np.tile(self.x_pos,96)
        self.ele_y = [i for i in range(6,3840, 20) for j in range(2)]
        self.ele_pos = np.column_stack((self.ele_x, self.ele_y))

# Read the bindary file "pre.bin"
def read_data(t0,tf):
    dt=1/2500.0
    m0=int(np.floor(t0/dt)-1)
    mf=int(np.ceil(tf/dt))
    mm = mf-m0
    n=384
    m = 2250000; 
    data=np.zeros((n,mm))
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/filtered_pre.bin", "rb") as fin:
        data[0,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='',offset=8*m0)
        for i in range(1,n):
            data[i,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='', offset=8*(m-mm))
    return data;


# Compute 2D kCSD in the entire space
def do_kcsd(pots,opts):
    length = pots.shape[1] # size of time dimensions
    pots = pots.reshape((len(opts.ele_pos), length)) 
    k = KCSD2D(opts.ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
               xmin=opts.xmin, xmax=opts.xmax,
               ymin=opts.ymin, ymax=opts.ymax,
               n_src_init=opts.n_src_init, src_type=opts.src_type,
               R_init=opts.R_init,lambd=opts.lambd,
               gdx=opts.gdx,gdy=opts.gdy)

    return k

# Compute 2D kCSD at the electrode positions and compute the error
# the opts need to contain the source locations:
# usually one should use opts.src_x=k.src_x; opts.src_x=k.src_y
def kcsd_error(pots,opts):
    ownk = oKCSD2D(opts.ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
                xmin=opts.xmin, xmax=opts.xmax,
                ymin=opts.ymin, ymax=opts.ymax,
                n_src_init=opts.n_src_init, src_type=opts.src_type,
                R_init=opts.R_init,lambd=opts.lambd,
                own_est=np.array((opts.ele_x,opts.ele_y)),own_src=(opts.src_x,opts.src_y)) 

    est_pots = ownk.values('POT')
    error = np.linalg.norm(ownk.pots-est_pots)**2 #+ opts.lambd*np.linalg.norm(est_pots)**2;
    return ownk,error;


# CV to obtain reasonable values of R and lambda.
# A bit heavy for long time series
def validate():
    opts=kcsd_opts()
    pots =read_data(100.0,100.0)
    lambdas=np.logspace(-12, -1, num=12)
    Rs=np.linspace(10, 100, num=10)
    hs=np.arange(1,101,10)
    cverrors=np.zeros([lambdas.size,Rs.size,hs.size])
    for i,h in enumerate(hs):
        opts.h=h
        k = do_kcsd(pots,opts)
        k.cross_validate(lambdas=lambdas, Rs=Rs)
        cverrors[:,:,i]=k.errs
    np.savetxt("/home/pclusella/Documents/Data/UPO-tACs/7_results/kCSD_validation/cv_errors.dat", 
               cverrors)
    return cverrors
    
def export_at_electrodes(opts):
    pots =read_data(0.0004,900.0)
    redk = oKCSD2D(opts.ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
            xmin=opts.xmin, xmax=opts.xmax,
            ymin=opts.ymin, ymax=opts.ymax,
            n_src_init=opts.n_src_init, src_type=opts.src_type, 
            R_init=opts.R_init,lambd=opts.lambd,
            own_est=np.array((opts.ele_x,opts.ele_y)),own_src=(opts.src_x,opts.src_y)) 
    est_csd = redk.values('CSD')
    # est_pot = redk.values('POT') # can export estimated potentials
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/kCSD_electrodes_pre.bin", "wb") as fout:
        est_csd.tofile(fout,"");
    return est_csd;

def export_at_centers(opts):
    pots =read_data(0.0004,900.0)
    loc_y=np.arange(opts.ymin,opts.ymax,opts.gdy)
    loc_x=np.zeros(loc_y.size)
    redk = oKCSD2D(opts.ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
            xmin=opts.xmin, xmax=opts.xmax,
            ymin=opts.ymin, ymax=opts.ymax,
            n_src_init=opts.n_src_init, src_type=opts.src_type, 
            R_init=opts.R_init,lambd=opts.lambd,own_est=np.array((loc_x,loc_y)),own_src=(opts.src_x,opts.src_y)) 
    est_csd = redk.values('CSD')
    # est_pot = redk.values('POT') # can export estimated potentials
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_Raw/kCSD_centers_pre.bin", "wb") as fout:
        est_csd.tofile(fout,"");
    return est_csd;



# validate()

# Requires a good amount of free RAM memory (tested in a system with 64Gb)
opts=kcsd_opts()
pots =read_data(100.0,100.0004)
k = do_kcsd(pots,opts)
opts.src_x = k.src_x
opts.src_y = k.src_y # this should be made differently

csd_electro=export_at_electrodes(opts)
csd_centers=export_at_centers(opts)


pots =read_data(100.0,100.0004)
k = do_kcsd(pots,opts)
# est_csd = k.values('CSD')


# -- 
# # Read s seconds, compute, and plot the kCSD
# s = 1.0
# pots =read_data(0.0004,0.0004)
# k = do_kcsd(pots,opts)
# est_csd = k.values('CSD')
# opts.src_x=k.src_x;
# opts.src_y=k.src_y;
# export_at_electrodes(opts)
# redk, err = kcsd_error(ele_pos,pots,ele_x,ele_y,opts)
# err
# # --

# # --
# est_csd = k.values('CSD')
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

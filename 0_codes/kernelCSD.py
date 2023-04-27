# Execute this with:
# exec(open("../0_codes/kernelCSD.py").read())
# opts=kcsd_opts()
# csd_centers=export_at_centers(opts,"pre")
import sys
import numpy as np
import matplotlib.pyplot as plt
from kcsd import KCSD2D
from kcsd import oKCSD2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Optimal paramters for cortex, e.g., 
# h=100, R=60 and l=1e-5
# h=10, R=60 and l=1e-6
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
def read_data(t0,tf,fl,m):
    dt = 1/2500.0
    m0 = int(np.floor(t0/dt)-1)
    mf = int(np.ceil(tf/dt))
    mm = mf-m0
    n = 136
    # m = 2250000; 
    # m = 9000000; 
    data=np.zeros((n,mm))
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_data/cortex_"+fl+".bin", "rb") as fin:
        data[0,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='',offset=8*m0)
        for i in range(1,n):
            data[i,:]=np.fromfile(fin, dtype=np.double, count=mm, sep='', offset=8*(m-mm))
    return data;

# Remove broken electrodes
def remove_broken(pots,opts,broken):
    # rmv=np.array([56,135,191,198,325])
    pots=np.delete(pots,broken,0)
    opts.ele_pos=np.delete(opts.ele_pos,broken,0)
    return pots,opts;

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

def select_cortex(pots,opts,n0,nf):
    # pots=pots[n0:nf];
    opts.ele_pos=opts.ele_pos[n0:nf+1];
    opts.ele_y=opts.ele_y[n0:nf+1];
    opts.ele_x=opts.ele_x[n0:nf+1];
    opts.ymin=np.floor(n0/2)*20-400.0
    opts.ymax=3840.0; #np.floor(nf/2)*20+400.0
    return pots,opts;


def export_at_centers(opts,fl,m):
    pots =read_data(0.0004,m*0.0004,fl,m)
    n0=226;
    nf=361;
    pots,opts=select_cortex(pots,opts,n0,nf)
    rmv=np.array([325])-n0
    pots,opts = remove_broken(pots,opts,rmv)
    loc_y=np.arange(np.floor(n0/2)*20,np.ceil(nf/2)*20,10)
    loc_x=np.zeros(loc_y.size)
    
    step=20
    x = np.arange(opts.xmin,opts.xmax+step,step)
    y = np.arange(opts.ymin,opts.ymax+step,step)
    
    xx=np.tile(x,y.size)
    yy=y.repeat(x.size)
    # sources = np.column_stack((xx,yy))

    redk = oKCSD2D(opts.ele_pos, pots, h=opts.h, sigma=opts.sigma,                                                                                                                                                       
        xmin=opts.xmin, xmax=opts.xmax,
        ymin=opts.ymin, ymax=opts.ymax,
        n_src_init=opts.n_src_init, src_type=opts.src_type, 
        R_init=opts.R_init,lambd=opts.lambd,own_est=np.array((loc_x,loc_y)),own_src=(xx,yy)) 
    est_csd = redk.values('CSD')
    # est_pot = redk.values('POT') # can export estimated potentials
    with open("/home/pclusella/Documents/Data/UPO-tACs/1_data/kcsd_"+fl+".bin", "wb") as fout:
        est_csd.tofile(fout,"");
    return est_csd;


def process_data(fl,m):
    opts=kcsd_opts()
    csd_centers=export_at_centers(opts,fl,m)
    return ;


mpre = int(sys.argv[1])
mpost= int(sys.argv[2])
process_data("pre",mpre);
process_data("post",mpost);

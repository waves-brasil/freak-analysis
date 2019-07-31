import numpy as np
import pandas as pd
import xarray as xr
import progressbar


import freak_utils as frk

# @profile
# def run():

fname = './data/SOFS-2-2011-11-08T220019.bin.nc'
# read data into dataset (SOFS-2 times are all good)
ds = frk.read_nc(fname)

df = pd.DataFrame(index=ds.TIME, columns={'hs', 'hmax', 'fw', 'hm0', 'tp'})

bar = progressbar.ProgressBar(maxval=len(ds.TIME), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.SimpleProgress()])
bar.start()
for i,it in enumerate(ds.TIME):
    bar.update(i+1)
    
    # z-acceleration to z-displacement
    zac = - ds.Acceleration.sel(TIME=it, vector=2) # invert for SOFS-2 (MRU mounted upside down)
    zdisp = frk.int2byf(-zac, dt=0.2, cut_off=1/0.03) # vertical displacement (eta)
    # horizontal accelerations for sanity check
    xac = ds.Acceleration.sel(TIME=it, vector=0)
    yac = ds.Acceleration.sel(TIME=it, vector=1)
    
    # sanity check
    if not (np.any(np.isnan(xac)) or np.any(np.isnan(yac)) or np.any(np.isnan(zac)) or xac.sum()==0 or yac.sum()==0 or zac.sum()==0 ):
        # spectrum (mostly for sanity check with wave parameters)
        spec = frk.spec1(zdisp, nfft=512, fs=5)
        f, S = spec[:,0], spec[:,1]
        tp = 1/f[S==max(S)][0]
        hm0 = 4*np.sqrt( sum(S) * (f[1]-f[0]) )
    else:
        continue

    if hm0 < 20 and tp < 25:
        izd, iza = frk.zero_cross(zdisp)
        ca_d, cr_d, h_d, ca_a, cr_a, cr, h_a, t_a = frk.ind_params(izd, iza, zdisp)
        hs, hmax, Fw = frk.calc_hs_fw(h_a)
        
        df.loc[it.values][['hs','hmax','fw','hm0','tp']] = [hs, hmax, Fw, hm0, tp]
    
bar.finish()

df.dropna(inplace=True)

# if __name__ == '__main__':
#     run()

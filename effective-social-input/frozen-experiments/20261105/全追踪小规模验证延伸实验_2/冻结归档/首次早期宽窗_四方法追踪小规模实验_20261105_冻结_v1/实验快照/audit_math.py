"""Independent scalar KF and neural-field replay; no simulation imports."""
import bootstrap
import math
import numpy as np
from pathlib import Path
import json
CFG=json.loads((Path(__file__).parent/"experiment_config.json").read_text(encoding="utf-8"))
N,K,DT=CFG["n"],CFG["steps"],CFG["dt"]

def circular(x):
    return np.arctan2(np.sin(x), np.cos(x))


def topology():
    selected=np.zeros((12,12),bool)
    hidden=selected.copy()
    for rank in range(11):
        i=rank+1
        peers=[1+(rank+j)%11 for j in range(1,7)]
        selected[i,[0]+peers]=True
        hidden[i,[0]+peers[:3]]=True
    return selected,hidden


def scalar_filter(observations,available,selected):
    """Independent scalar covariance update, with both Hold and KF candidates."""
    states=np.zeros((K,N,N,2))
    cov=np.zeros((K,N,N,2,2))
    last_values=np.full((K,N,N),np.nan)
    for i,j in zip(*np.nonzero(selected)):
        a=w=p00=p01=p10=p11=0.
        initialized=False
        last=float('nan')
        for k in range(K):
            if initialized:
                a+=DT*w
                p00,p01,p10,p11=(p00+DT*(p01+p10)+DT*DT*p11+1e-4*DT**3/3,
                                  p01+DT*p11+1e-4*DT**2/2,
                                  p10+DT*p11+1e-4*DT**2/2,p11+1e-4*DT)
            if available[k,i,j]:
                y=float(observations[k,i,j])
                if initialized:
                    innovation=math.atan2(math.sin(y-a),math.cos(y-a))
                    k0,k1=p00/(p00+1e-4),p10/(p00+1e-4)
                    a,w=a+k0*innovation,w+k1*innovation
                    p00,p01,p10,p11=((1-k0)*p00,(1-k0)*p01,p10-k1*p00,p11-k1*p01)
                else:
                    a,w,p00,p01,p10,p11=y,0.,1e-4,0.,0.,1.
                    initialized=True
                last=y
            a%=2*math.pi
            states[k,i,j]=[a,w]
            cov[k,i,j]=[[p00,p01],[p10,p11]]
            last_values[k,i,j]=last
    return states,cov,last_values


def control_replay(run):
    grid=np.arange(100)*2*np.pi/100
    diff=circular(grid[:,None]-grid[None,:])
    ring=np.cos(np.pi*np.sqrt(np.abs(diff)/np.pi))
    u=run['initial_u'][:,1:].copy()
    maximum=0.
    for k in range(K):
        weights=.012*run['edge_scales'][1:]*run['active'][k,1:]
        differences=circular(run['alpha'][1:,None,:]-run['centers'][k,1:,:,None])
        field=np.sum(weights[:,:,None]*np.exp(-.5*(differences/.4)**2),axis=1).T
        u=u+.3*(-u+ring@np.tanh(1000*u)/100+field)
        activation=np.maximum(np.tanh(1000*u),0.)
        v=.05*np.c_[(activation*np.cos(run['alpha'][1:].T)).sum(0),
                     (activation*np.sin(run['alpha'][1:].T)).sum(0)]
        maximum=max(maximum,float(np.abs(v-run['velocities'][k,1:]).max()))
    maximum=max(maximum,float(np.abs(u-run['final_u'][:,1:]).max()))
    assert maximum<1e-8, ('neural field replay',maximum)
    return maximum


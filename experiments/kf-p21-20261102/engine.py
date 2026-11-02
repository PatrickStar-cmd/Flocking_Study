"""Instrumented, batched stock CV/GM recursion; no truth enters this module.

Covariance is shared across edges observing the same mask. Batch dimension is
independent candidate windows. Output includes angular velocity (passive only).
"""
import numpy as np


def matrices(dt,q,tau=None):
    if tau is None or tau==float('inf'):
        return np.array([[1.,dt],[0.,1.]]),q*np.array([[dt**3/3,dt**2/2],[dt**2/2,dt]])
    rho=np.exp(-dt/tau)
    F=np.array([[1.,tau*(1-rho)],[0.,rho]])
    Q=q*np.array([[tau**2*dt-tau**3*(2*(1-rho)-(1-rho**2)/2),tau**2*(1-rho)**2/2],
                  [tau**2*(1-rho)**2/2,tau*(1-rho**2)/2]])
    return F,(Q+Q.T)/2


def run(observations,missing,dt=.3,q=1e-4,R=1e-4,tau=None,initial=None,hold=False):
    y=np.asarray(observations,dtype=float)
    if y.ndim==2:y=y[:,None,:]
    T,B,E=y.shape
    assert T>0 and E>0 and np.isfinite(y).all() and dt>0 and R>0 and q>=0
    missing=np.asarray(missing,dtype=bool)
    assert missing.shape==(T,)
    F,Q=matrices(dt,q,tau)
    theta=np.empty((T,B,E));omega=np.empty_like(theta);nis=np.full_like(theta,np.nan)
    P_record=np.empty((T,B,2,2));age=np.empty(T)
    if initial is None:
        assert not missing[0]
        x=np.stack([np.mod(y[0],2*np.pi),np.zeros((B,E))],axis=-1)
        P=np.broadcast_to(np.diag([R,1.]),(B,2,2)).copy()
        start=0
    else:
        x=np.stack([initial['theta'],initial['omega']],axis=-1).copy()
        P=initial['P'].copy(); start=-1
    a=0.
    for t in range(T):
        if hold:
            if not missing[t]:x[...,0]=np.mod(y[t],2*np.pi)
            a=a+dt if missing[t] else 0.
        elif t!=start:
            x=x@F.T
            P=np.matmul(F,np.matmul(P,F.T))+Q
            a+=dt
            if not missing[t]:
                innovation=np.arctan2(np.sin(y[t]-x[...,0]),np.cos(y[t]-x[...,0]))
                S=P[:,0,0]+R
                K=P[:,:,0]/S[:,None]
                nis[t]=innovation**2/S[:,None]
                x+=innovation[...,None]*K[:,None,:]
                A=np.broadcast_to(np.eye(2),(B,2,2)).copy()
                A[:,:,0]-=K
                P=np.matmul(A,np.matmul(P,A.swapaxes(-1,-2)))+R*K[:,:,None]*K[:,None,:]
                a=0.
            x[...,0]=np.mod(x[...,0],2*np.pi)
            P=(P+P.swapaxes(-1,-2))/2
        theta[t]=x[...,0];omega[t]=x[...,1]
        P_record[t]=P;age[t]=a
    if hold:
        omega[:]=np.nan;P_record[:]=np.nan
    return {'theta':theta,'omega':omega,'P':P_record,'nis':nis,'age':age}


def archive_branch(path,branch,fresh,indices,starts,initial,**metadata):
    """Exact float-bit XOR compression against shared all-fresh records.

    No information loss: archived uint64 xor shared uint64 recovers the exact
    float64 output. Pre200 is in shared fresh records; recovery400 is explicit.
    """
    assert not path.exists()
    arrays={}
    for name in ['theta','omega','nis','P']:
        baseline=fresh[name][indices,0]
        current=branch[name]
        assert baseline.shape==current.shape
        arrays[name+'_xor']=np.bitwise_xor(np.ascontiguousarray(current).view(np.uint64),np.ascontiguousarray(baseline).view(np.uint64))
    np.savez_compressed(path,**arrays,age=branch['age'],starts=starts,
        initial_theta=initial['theta'],initial_omega=initial['omega'],initial_P=initial['P'],
        **metadata)


def recover_branch(path,fresh):
    with np.load(path,allow_pickle=False) as data:
        starts=data['starts'].astype(int)
        K=data['theta_xor'].shape[0]
        indices=np.arange(K)[:,None]+starts[None,:]-1
        result={name:np.bitwise_xor(data[name+'_xor'],np.ascontiguousarray(fresh[name][indices,0]).view(np.uint64)).view(np.float64)
                for name in ['theta','omega','nis','P']}
        result['age']=data['age']
        return result

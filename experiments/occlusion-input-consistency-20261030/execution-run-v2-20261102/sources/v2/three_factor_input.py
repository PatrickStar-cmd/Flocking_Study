"""Early-v3 positive angular input: strength x count x relative-bearing weight."""
import numpy as np

def social_input(weights,bearings,heading,kappa=1.0):
    w=np.asarray(weights,float);theta=np.asarray(bearings,float);phi=np.asarray(heading,float)
    assert w.ndim==2 and theta.shape==w.shape and phi.shape==(w.shape[0],)
    assert np.isfinite(kappa) and kappa>0
    assert np.isfinite(w).all() and (w>=0).all() and np.isfinite(phi).all()
    active=w>0;assert np.isfinite(theta[active]).all()
    delta=np.where(active,theta-phi[:,None],0)
    q=np.exp(kappa*(np.cos(delta)-1));contribution=w*q
    A=w.sum(axis=1);S=contribution.sum(axis=1);n=active.sum(axis=1)
    C=np.divide(S,A,out=np.zeros_like(S),where=A>0)
    mean_weight=np.divide(A,n,out=np.zeros_like(A),where=n>0)
    return dict(S=S,A=A,n=n,mean_weight=mean_weight,angle_weight=C,
        Vx=np.sum(contribution*np.cos(delta),axis=1),Vy=np.sum(contribution*np.sin(delta),axis=1))

def checkpoint_input(checkpoint,remaining,model,kappa=1.0):
    """Use only the visible subset, with the common event heading as reference."""
    n=int(model['N']);k=int(model['full_neighbors']);r=int(remaining)
    assert 0<=r<=k
    pos=np.asarray(checkpoint['position'],float);head=np.asarray(checkpoint['heading'],float).reshape(n)
    order=np.asarray(checkpoint['maskOrder'],int)-1
    assert order.shape==(n,k) and all(set(row)==set(range(k)) for row in order)
    ids=(np.arange(n)[:,None]+np.arange(1,k+1)[None,:])%n
    dx=pos[0,ids]-pos[0,:,None];dy=pos[1,ids]-pos[1,:,None];L=model['L']
    dx-=L*np.round(dx/L);dy-=L*np.round(dy/L)
    visible=np.zeros((n,k),bool)
    for i in range(n):visible[i,order[i,:r]]=True
    observed=np.where(visible,np.arctan2(dy,dx),np.nan)
    weights=model['single_edge_weight']*visible
    return social_input(weights,observed,head,kappa)

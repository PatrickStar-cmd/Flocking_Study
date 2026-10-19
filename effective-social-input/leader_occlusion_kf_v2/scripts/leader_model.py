"""Protocol-faithful simulation, extending the SHA-verified reference implementation.

Control receives only sanitized bearings. Geometry and task evaluation are kept
in the simulation/sensor layer. Full state and selected-edge traces are saved.
"""
from dataclasses import asdict
import copy
import numpy as np
import reference_run_experiment as ref

MODES=['fresh','drop','drop-normalized','hold','linear-trend','kf-mean','kf-uncertainty']
DEVELOPMENT=[2026101401,2026101402,2026101403]
EVALUATION=list(range(2026101501,2026101531))

def config(n=20):
    return ref.Config(n=n,k=15,lead_steps=100,total_steps=800,cue_amplitude=.32,
                      cue_shape='cosine',geometry_scale=1)

def leader_count(c): return 3 if c.n==20 else int(.15*c.n)

def topology(c):
    if c.n==20:return ref.topology(c)
    n=c.n; L=leader_count(c); F=n-L; s=np.zeros((n,n),bool)
    for f in range(F):
        s[L+f,(f+np.arange(3))%L]=True
        s[L+f,L+(f+np.arange(1,13))%F]=True
    for l in range(L):
        s[l,(l+np.arange(1,3))%L]=True
        s[l,L+(13*l+np.arange(13))%F]=True
    h=s.copy();h[:L]=False;h[:,L:]=False
    assert np.all(s.sum(1)==15) and not np.diag(s).any() and h.sum()==3*F
    return s,h

def initialize(seed,c):
    if c.n==20:return ref.initialize(seed,c)
    L=leader_count(c);F=c.n-L;r=int(np.ceil(np.sqrt(F)))
    rng=np.random.Generator(np.random.PCG64(seed));xy=np.empty((c.n,2))
    xy[:L,0]=492;xy[:L,1]=500+8*(np.arange(L)-(L-1)/2)
    j=np.arange(F);xy[L:,0]=480-8*(j//r);xy[L:,1]=500+8*((j%r)-(r-1)/2)
    xy+=rng.normal(0,.4,(c.n,2));heading=rng.normal(0,.03,c.n)
    alpha=(np.arange(c.m)[None,:]*2*np.pi/c.m+heading[:,None])%(2*np.pi)
    u=.1*np.cos(alpha.T-heading[None,:])+rng.normal(0,.001,(c.m,c.n))
    noise=rng.normal(0,c.noise_std,(c.total_steps,c.n,c.n))
    return xy,heading,alpha,u,noise

class ForecastMemory(ref.Memory):
    """Shadow KF unchanged; linear baseline uses only contiguous legal history."""
    def __init__(self,selected,cfg):
        super().__init__(selected,cfg)
        self.history=np.full((10,*selected.shape),np.nan)
        self.count=np.zeros(selected.shape,int)
        self.previous_fresh=np.zeros_like(selected)
        self.slope=np.zeros(selected.shape)

    def input(self,observed,fresh,alpha,heading,mode):
        if mode!='linear-trend':return super().input(observed,fresh,alpha,heading,mode)
        restarted=fresh & ~self.previous_fresh
        self.history[:,restarted]=np.nan;self.count[restarted]=0
        self.history[1:,fresh]=self.history[:-1,fresh]
        self.history[0,fresh]=observed[fresh]
        self.count[fresh]=np.minimum(self.count[fresh]+1,10)
        # Unwrap backwards from the newest observation. No predicted value enters history.
        h=self.history;dif=ref.wrap(h[:-1]-h[1:])
        y=np.concatenate([h[:1],h[:1]-np.cumsum(dif,axis=0)],axis=0)
        valid=np.arange(10)[:,None,None]<self.count[None]
        y=np.where(valid,y,0.);t=-np.arange(10)[:,None,None]*self.cfg.dt
        cnt=np.maximum(self.count,1);st=np.sum(t*valid,axis=0);sy=np.sum(y,axis=0)
        numerator=np.sum(t*y,axis=0)-st*sy/cnt
        denominator=np.sum(t*t*valid,axis=0)-st*st/cnt
        self.slope=np.divide(numerator,denominator,out=np.zeros_like(numerator),where=denominator>0)
        self.previous_fresh=fresh.copy()
        _,d=super().input(observed,fresh,alpha,heading,'hold')
        active=d['active'];safe_age=np.where(self.initialized,self.age,0.)
        center=np.where(fresh,observed,(self.last+safe_age*self.slope)%(2*np.pi))
        center=np.where(active,center,0.);weights=d['weights']
        kernels=np.exp(-.5*(ref.wrap(alpha[:,None,:]-center[:,:,None])/self.cfg.sigma)**2)
        field=np.sum(weights[:,:,None]*kernels,axis=1).T
        factor=np.exp(np.cos(center-heading[:,None])-1.)
        S=np.sum(weights*factor,axis=1);A=d['A']
        d.update(center=center,S=S,C=np.divide(S,A,out=np.zeros_like(S),where=A>0),
                 S_memory=np.sum(weights*factor*(active & ~fresh),axis=1))
        return field,d

def snapshot(memory,xy,heading,alpha,u):
    return dict(xy=xy.copy(),heading=heading.copy(),alpha=alpha.copy(),u=u.copy(),
        filter_x=memory.x.copy(),filter_P=memory.P.copy(),initialized=memory.initialized.copy(),
        last_bearing=memory.last.copy(),age=memory.age.copy(),history=memory.history.copy(),
        history_count=memory.count.copy(),previous_fresh=memory.previous_fresh.copy(),slope=memory.slope.copy())

def simulate(seed,c,mode,duration=60,start=200,initial=None):
    initial=initialize(seed,c) if initial is None else initial
    xy,heading,alpha,u,noise=[a.copy() for a in initial]
    selected,hidden=topology(c);L=leader_count(c);ii,jj=np.where(selected);ei=len(ii)
    memory=ForecastMemory(selected,c);W=ref.ring(c.m);rows=[];states={}
    states['xy']=np.empty((801,c.n,2));states['headings']=np.empty((801,c.n))
    states['u']=np.empty((801,c.m,c.n));states['velocity']=np.empty((800,c.n,2))
    states['xy'][0]=xy;states['headings'][0]=heading;states['u'][0]=u
    names=['observed','true','center','weight','filter_theta','filter_omega','P11','P12','P22','age','active','fresh']
    trace={k:np.empty((800,ei),dtype=bool if k in ['active','fresh'] else float) for k in names}
    agent={k:np.empty((800,c.n)) for k in ['S','S_oracle','A','C']}
    states['social_field']=np.empty((800,c.m,c.n))
    fixed_nearest=None;checkpoint={};lasttruth=None
    for step in range(c.total_steps):
        theta,distance=ref.geometry(xy,c)
        if step==start:
            checkpoint=snapshot(memory,xy,heading,alpha,u)
            d0=distance.copy();np.fill_diagonal(d0,np.inf);fixed_nearest=np.argsort(d0,axis=1)[:,:5]
        blackout=mode!='fresh' and start<=step<start+duration
        fresh=selected & ~(hidden if blackout else np.zeros_like(hidden))
        observed=np.where(fresh,(theta+noise[step])%(2*np.pi),np.nan)
        field,d=memory.input(observed,fresh,alpha,heading,mode)
        states['social_field'][step]=field
        command=min(c.turn_angle,c.turn_rate*max(0,(step-c.warm_steps+1)*c.dt))
        field[:,:L]+=c.cue_amplitude*np.cos(alpha[:L]-command).T
        oracle=c.h0*np.sum(selected*np.exp(np.cos(theta-heading[:,None])-1.),axis=1)
        for key,val in [('observed',observed),('true',theta),('center',d['center']),('weight',d['weights']),
            ('filter_theta',memory.x[:,:,0]),('filter_omega',memory.x[:,:,1]),('P11',memory.P[:,:,0,0]),
            ('P12',memory.P[:,:,0,1]),('P22',memory.P[:,:,1,1]),('age',memory.age),('active',d['active']),('fresh',fresh)]:trace[key][step]=val[ii,jj]
        for key in ['S','A','C']:agent[key][step]=d[key]
        agent['S_oracle'][step]=oracle
        prediction_mask=hidden & blackout;active_mask=prediction_mask & d['active']
        shadow=ref.wrap(memory.x[:,:,0]-theta)[prediction_mask];V=memory.P[:,:,0,0][prediction_mask]
        u+=c.dt*(-u+W@np.tanh(c.beta*u)/c.m+field)
        positive=np.maximum(np.tanh(c.beta*u),0.)
        v=c.speed_scale*np.c_[np.sum(positive*np.cos(alpha.T),axis=0),np.sum(positive*np.sin(alpha.T),axis=0)]
        speed=np.linalg.norm(v,axis=1);moving=speed>1e-12
        direction=np.divide(v,speed[:,None],out=np.zeros_like(v),where=moving[:,None])
        heading[moving]=np.arctan2(v[moving,1],v[moving,0])%(2*np.pi)
        xy=(xy+c.dt*v)%c.arena
        states['xy'][step+1]=xy;states['headings'][step+1]=heading
        states['u'][step+1]=u;states['velocity'][step]=v
        _,after=ref.geometry(xy,c);nearest=np.argsort(after,axis=1)[:,:6]
        noself=after.copy();np.fill_diagonal(noself,np.inf)
        leader_z=np.mean(np.exp(1j*heading[:L]));follower_z=np.mean(np.exp(1j*heading[L:]))
        e=np.abs(ref.wrap(heading-command));err=ref.wrap(d['center']-theta)
        gated_age=prediction_mask & memory.initialized & (memory.age>c.max_age)
        gated_variance=prediction_mask & memory.initialized & (memory.P[:,:,0,0]>c.max_variance)
        row=dict(step=step,time=(step+1-start)*c.dt,command_deg=np.degrees(command),blackout=int(blackout),
            GO=np.linalg.norm(direction.sum(0))/c.n,GO_followers=np.linalg.norm(direction[L:].sum(0))/(c.n-L),
            GO_leaders=np.linalg.norm(direction[:L].sum(0))/L,
            LO=np.linalg.norm(direction[nearest].sum(axis=1),axis=1).mean()/6,
            fixed5_distance=after[np.arange(c.n)[:,None],fixed_nearest].mean() if fixed_nearest is not None else np.nan,
            nearest_distance=noself.min(axis=1).mean(),pair_distance=after[np.triu_indices(c.n,1)].mean(),
            leader_follower_distance=after[L:,:L].mean(),minimum_distance=noself.min(),
            follower_turn_error_deg=np.degrees(e[L:]).mean(),leader_turn_error_deg=np.degrees(e[:L]).mean(),
            follower_leader_error_deg=np.degrees(np.abs(ref.wrap(heading[L:]-np.angle(leader_z)))).mean(),
            follower_turn_angle_deg=np.degrees(np.angle(follower_z)),leader_turn_angle_deg=np.degrees(np.angle(leader_z)),
            follower_resultant=abs(follower_z),leader_resultant=abs(leader_z),mean_speed=speed.mean(),moving_fraction=moving.mean(),
            S=d['S'][L:].mean(),A=d['A'][L:].mean(),C=d['C'][L:].mean(),S_memory=d['S_memory'][L:].mean(),
            S_oracle=oracle[L:].mean(),S_abs_error=np.abs(d['S'][L:]-oracle[L:]).mean(),
            S_hidden_oracle=(c.h0*np.sum(hidden*np.exp(np.cos(theta-states['headings'][step,:,None])-1.),axis=1))[L:].mean(),
            fresh_edges=int(fresh.sum()),active_edges=int(d['active'].sum()),hidden_active_edges=int(active_mask.sum()),
            gated_age_edges=int(gated_age.sum()) if mode in ['hold','linear-trend','kf-mean','kf-uncertainty'] else 0,
            gated_variance_edges=int(gated_variance.sum()) if mode.startswith('kf') else 0,
            bearing_error_deg=np.degrees(np.abs(err[active_mask])).mean() if active_mask.any() else np.nan,
            shadow_kf_rmse_rad=np.sqrt(np.mean(shadow**2)) if shadow.size else np.nan,
            shadow_kf_coverage=np.mean(np.abs(shadow)<=1.96*np.sqrt(V)) if V.size else np.nan,
            hidden_bearing_change_deg=np.degrees(np.abs(ref.wrap(theta-lasttruth))[hidden]).mean() if step>=start else np.nan)
        if step==start-1:lasttruth=theta.copy()
        rows.append(row)
    return dict(rows=rows,states=states,edges=trace,agent=agent,alpha=alpha,selected=selected,hidden=hidden,
        edge_i=ii,edge_j=jj,fixed_nearest=fixed_nearest,checkpoint=checkpoint,
        final=snapshot(memory,xy,heading,alpha,u))

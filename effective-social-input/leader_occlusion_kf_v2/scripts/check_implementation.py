"""Meaningful reference, interface, estimator and protocol checks before running."""
from pathlib import Path
import json,time
import numpy as np
from threadpoolctl import threadpool_limits
import leader_model as model
import reference_run_experiment as ref

ROOT=Path(__file__).resolve().parent

def run_checks(out):
    out.mkdir(parents=True,exist_ok=True);records=[]
    def check(name,fn):
        start=time.perf_counter();detail=fn();records.append(dict(check=name,status='PASS',seconds=time.perf_counter()-start,detail=detail))
        print('PASS',name,flush=True)
    c=model.config();selected,hidden=model.topology(c)
    def graphs():
        from scipy.sparse.csgraph import connected_components
        result=[]
        for n in [20,40,60]:
            cfg=model.config(n);s,h=model.topology(cfg);L=model.leader_count(cfg)
            assert s.sum()==15*n and h.sum()==3*(n-L)
            assert np.all(s[L:,:L].sum(1)==3)
            assert connected_components(s,directed=True,connection='strong')[0]==1
            result.append(dict(N=n,leaders=L,edges=int(s.sum()),hidden=int(h.sum())))
        return result
    check('graph_counts_and_strong_connectivity',graphs)
    def scalar():
        cfg=model.config();s=np.zeros((20,20),bool);s[0,1]=True;mem=model.ForecastMemory(s,cfg)
        F=np.array([[1,.3],[0,1]]);Q=cfg.q*np.array([[.3**3/3,.3**2/2],[.3**2/2,.3]])
        x=None;P=None;mx=0.
        for k in range(80):
            z=(6.1+.025*k)%(2*np.pi);fresh=s.copy();obs=np.full((20,20),np.nan)
            if 30<=k<55:fresh[:]=False
            else:obs[0,1]=z
            if x is None:x=np.array([z,0.]);P=np.diag([cfg.R,1.])
            else:
                x=F@x;P=F@P@F.T+Q
                if fresh[0,1]:
                    K=P[:,0]/(P[0,0]+cfg.R);x+=K*ref.wrap(z-x[0]);J=np.eye(2)-np.outer(K,[1,0])
                    P=J@P@J.T+cfg.R*np.outer(K,K)
                x[0]%=2*np.pi
            mem.input(obs,fresh,np.tile(np.arange(100)*2*np.pi/100,(20,1)),np.zeros(20),'kf-mean')
            diff=max(np.max(np.abs(x-mem.x[0,1])),np.max(np.abs(P-mem.P[0,1])))
            mx=max(mx,float(diff));assert diff<1e-12
        assert np.linalg.eigvalsh(mem.P[0,1]).min()>0
        return dict(maxDifference=mx)
    check('scalar_KF_wrap_missing_Joseph',scalar)
    def interface():
        cfg=model.config();s,h=model.topology(cfg);initial=model.initialize(2026101401,cfg);alpha=initial[2];heading=initial[1]
        mem=model.ForecastMemory(s,cfg);obs=np.where(s,.4,np.nan);mem.input(obs,s,alpha,heading,'kf-mean')
        import copy
        fresh=s & ~h;observed=np.where(fresh,.5,np.nan)
        a=copy.deepcopy(mem);b=copy.deepcopy(mem)
        # Evaluation truths are absent from this controller interface.
        truth1=np.full((20,20),.3);truth2=np.full((20,20),2.9)
        fa,da=a.input(observed,fresh,alpha,heading,'kf-mean');fb,db=b.input(observed,fresh,alpha,heading,'kf-mean')
        assert np.array_equal(fa,fb) and np.array_equal(a.x,b.x)
        bad=observed.copy();bad[h]=truth2[h]
        try:mem.input(bad,fresh,alpha,heading,'kf-mean')
        except AssertionError:pass
        else:raise AssertionError('finite hidden observation accepted')
        cold=model.ForecastMemory(s,cfg);f,d=cold.input(np.full((20,20),np.nan),np.zeros_like(s),alpha,heading,'kf-mean')
        assert not d['active'].any() and not f.any()
        for mode,attr,value in [('hold','age',121.),('kf-mean','variance',2.)]:
            m=copy.deepcopy(a)
            if attr=='age':m.age[h]=value
            else:m.P[:,:,0,0][h]=value
            _,d=m.input(observed,fresh,alpha,heading,mode);assert not d['active'][h].any()
            _,d=m.input(obs,s,alpha,heading,mode);assert d['active'][s].all() and np.all(m.age[s]==0)
        return 'Hidden truth interface, cold start, age/variance gate and recovery checked'
    check('causal_interface_and_gates',interface)
    def trend():
        cfg=model.config();s=np.zeros((20,20),bool);s[0,1]=True;alpha=np.tile(np.arange(100)*2*np.pi/100,(20,1))
        for omega in [0,.03]:
            mem=model.ForecastMemory(s,cfg)
            for k in range(25):
                obs=np.full((20,20),np.nan);fresh=s.copy()
                if k<15:obs[0,1]=(6.15+omega*k*.3)%(2*np.pi)
                else:fresh[:]=False
                _,d=mem.input(obs,fresh,alpha,np.zeros(20),'linear-trend')
                if k>=15:assert abs(ref.wrap(d['center'][0,1]-(6.15+omega*k*.3)))<1e-12
        return 'Zero and constant rate, wrap crossing, no forecast writes to history'
    check('linear_trend_analytic',trend)
    def uncertainty():
        from numpy.polynomial.hermite import hermgauss
        cfg=model.config();s=np.zeros((20,20),bool);s[0,1]=True
        alpha=np.tile(np.arange(100)*2*np.pi/100,(20,1));mem=model.ForecastMemory(s,cfg)
        obs=np.full((20,20),np.nan);obs[0,1]=.7;mem.input(obs,s,alpha,np.zeros(20),'kf-uncertainty')
        mem.P[0,1]=np.array([[.09,0],[0,0]]);mem.x[0,1,1]=0
        f,d=mem.input(np.full((20,20),np.nan),np.zeros_like(s),alpha,np.zeros(20),'kf-uncertainty')
        x,w=hermgauss(96);V=mem.P[0,1,0,0];mu=mem.x[0,1,0]
        numeric=np.sum(w*np.exp(np.cos(mu+np.sqrt(2*V)*x)-1))/np.sqrt(np.pi)
        assert abs(d['S'][0]/cfg.h0-numeric)<1e-12
        g=np.exp(-.5*(ref.wrap(alpha[0]-mu)/.4)**2)
        assert abs(f[:,0].sum()/cfg.h0-g.sum())<1e-10 and f.min()>=0
        return dict(SQuadratureDifference=float(abs(d['S'][0]/cfg.h0-numeric)))
    check('uncertainty_S_independent_quadrature_and_mass',uncertainty)
    def reference_parity():
        cfg=model.config();maxima={};initial=model.initialize(2026101401,cfg)
        for mode in ref.MODES:
            a=model.simulate(2026101401,cfg,mode,60,initial=initial)
            b=ref.simulate(2026101401,cfg,mode,60,initial=initial)
            err=max(np.max(np.abs(a['states']['xy']-b['xy'])),np.max(np.abs(a['states']['u'][-1]-b['u'])))
            assert err<1e-10,(mode,err)
            for metric in ['GO','LO','fixed5_distance','follower_turn_error_deg','S','S_abs_error']:
                av=np.array([r[metric] for r in a['rows']]);bv=np.array([r[metric] for r in b['rows']])
                assert np.allclose(av,bv,rtol=0,atol=1e-10,equal_nan=True),(mode,metric)
            maxima[mode]=float(err)
        return maxima
    check('six_methods_protocol_reference_parity',reference_parity)
    def no_occlusion():
        cfg=model.config();initial=model.initialize(2026101402,cfg)
        fresh=model.simulate(2026101402,cfg,'fresh',0,initial=initial)
        for mode in model.MODES[1:]:
            r=model.simulate(2026101402,cfg,mode,0,initial=initial)
            for field in ['xy','headings','u','velocity','social_field']:assert np.array_equal(fresh['states'][field],r['states'][field]),(mode,field)
        return 'All seven methods exactly equal when T=0'
    check('T0_seven_methods',no_occlusion)
    (out/'software_checks.json').write_text(json.dumps(dict(status='PASS',checks=records),indent=2),encoding='utf-8')
    return records

if __name__=='__main__':
    with threadpool_limits(limits=1):run_checks(ROOT/'checks/preflight')

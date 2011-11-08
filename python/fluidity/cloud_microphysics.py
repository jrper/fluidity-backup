#default_physical_parameters={'m':0.018}
import numpy as np
import fluidity_tools

class microphysics_forcing(dict):
    
    def __init__(self,forcing,args,bounds=()):
        dict.__init__(self)
        self.forcing=forcing
        self.bounds=bounds
        for arg in args:
            try:
                self[arg[0]]=arg[1]
            except TypeError:
                self[arg]=1.0

    def bound(self,f,alpha):
        if not self.bounds==():
            lbound=self.bounds[0]
            ubound=self.bounds[1]
            return np.where(f>0.0,f*alpha[ubound],f*alpha[lbound])
        else:
            return f

class microphysics_model(object):

    """ An example factory class for python microphysics models. The state passed as input has the source terms of the forced fields set to match the forcings provided, subject to imposed constraints on the water phases, namely that they lie between 0 and 1."""

    def __init__(self,state,
                 prescribed_fields={},
                 forced_fields=tuple(),
                 python_fields=tuple(),
                 forcings=tuple(),
                 slip_vels=tuple(),
                 dt=None):
        self.state=state
#        print state.scalar_fields
        self.prescribed_fields=prescribed_fields
        self.forced_fields=forced_fields.keys()
        self.dt=dt
        self.fields={}
        self.oldfields={}
        self._get_fields(prescribed_fields,forced_fields)
        self._add_python_fields(python_fields)
        self._calculateForcings()
        self._calculateSlipVels()
    
    def _get_fields(self,prescribed_fields,forced_fields):
        for field,name in prescribed_fields.items():
            self.fields[field]=self.state.scalar_fields[name].val[:]
            self.oldfields[field]=self.state.scalar_fields['Old'+name].val[:]
        for field,name in forced_fields.items():
            self.fields[field]=self.state.scalar_fields['Iterated'+name].val[:]
            self.oldfields[field]=self.state.scalar_fields['Old'+name].val[:]
        for field,name in forced_fields.items():
            self.fields['delta'+field]=(
                self.state.scalar_fields[name+'Source'].val)
            self.fields['delta'+field][:]=0
            try:
                self.fields[field+'slip']=(
                    self.state.scalar_fields[name+'SinkingVelocity'].val)
            except:
                pass
        self.fields['dt']=self.dt
        self.oldfields['dt']=self.dt


    def _add_python_fields(self,python_fields):
        for field_name,func in python_fields:
            self.fields[field_name]=func(**self.fields)
            self.oldfields[field_name]=func(**self.oldfields)

    def _calculateForcings(self):

        a=0.0

        alpha={}
        dF={}

        for F in self.forcings:
#            f=F.forcing(**self.fields)
            oldf=F.forcing(**self.oldfields)
            for field_name,scale in F.items():
                if field_name in dF:
                    dF[field_name]+=scale*(oldf)
                else:
                    dF[field_name]=scale*(oldf)

        print dF.keys()

        for field in dF:
            fv=self.fields[field]
            dF[field]=np.where(abs(dF[field])<1.0e-8,fv/self.dt,dF[field]) 
            A=-fv/(self.dt*dF[field])
            A=np.where(np.isnan(A),1.0,A)
            alpha[field]=np.where((A>1.0)*(fv>0.0),1.0,A)
            alpha[field]=np.where((A<0.0)*(fv>0.0),1.0,alpha[field])
            alpha[field]=np.where((A>0.0)*(fv<0.0),1.0,alpha[field])
            alpha[field]=np.where((A<0.0)*(fv<0.0),
                                  np.minimum(-A,1.0e6),alpha[field])

        for F in self.forcings:
#            f=F.bound(F.forcing(**self.fields),alpha)
            oldf=F.bound(F.forcing(**self.oldfields),alpha)
            for field_name,scale in F.items():
                q_w=self.fields['q_w']
                if field_name[0]=='q':
                    self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+fix(scale*(oldf),-q_w,q_w))
                else:
                    self.fields['delta'+field_name][:]=(self.fields['delta'+field_name][:]+scale*(oldf))   

    def _calculateSlipVels(self):
        pass


L_v=2500000.0
#m=0.018
m=3.0e-3
rho_w=1020.0
R_v=461.5
D_v=2.0e-5
K_a=0.03
N_0r=10.0e7
a_r=201
nu=20.0
rho_0=1.0
k_c=1e-3
a_c=5e-4

def fix(x,a,b):
    return np.maximum(np.minimum(x,b),a)

def testing(state,dt):
    M=hot_moist_microphysics_model(state,dt)


        


        
def get_gas_density(rho=None,q_g=None,q_v=None,q_c=None,q_r=None,**kwargs):
#    print 'Gas Density'
    rho_g=rho*q_g/(1.0-rho*(q_c+q_r)/rho_w)
    return np.where(rho_g>0,rho_g,np.inf)

def get_q_g(q_v,q_c,q_r,**kwargs):
    return 1.0-q_c-q_r

def get_q_w(q_v,q_c,q_r,**kwargs):
    return q_v+q_c-q_r

def get_q_sat(rho=None,p=None,T=None,rho_g=None,q_v=None,
              q_c=None,q_r=None,**kwargs):
    q_g=1.0-q_r-q_c
    print 'T : ', np.min(T), np.max(T)
    T=np.where(T<273.0,273.0,T)
    P_sat=611.2*np.exp(17.62*(T-273.0)/(T-30.0))
    q_sat =P_sat/R_v*(q_g/(rho_g*T))

    return q_sat
    
def supersaturation_fix(q_v=None,q_sat=None,dt=None,**kwargs):
    return np.maximum((q_v-q_sat)/dt,0.0)


def cloudwater_condensation(rho=None,p=None,T=None,
                            q_v=None,q_c=None,q_r=0.0,
                            q_sat=None,rho_g=None,**kwargs):

    from numpy import pi

    S_w=np.minimum(q_v/q_sat,1.0)
    q_g=1.0-q_c-q_r
    r_c=10.0e-6
    N_c=3.0*q_c/(4.0*rho_w*pi*r_c**3)
#    N_c=10000

    F_k=(L_v/(R_v*T)-1.0)*(L_v/(K_a*T))

#    P_sat=np.exp(20.386-5132.0/T)

    P_sat=611.2*np.exp(17.62*(T-273.0)/(T-30.0))

    F_d=(R_v*T)/(P_sat*D_v)

    ConC=4.0*np.pi*(q_g/rho_g)*N_c*r_c*(S_w-1.0)/(F_k+F_d)
#    ConC=4.0*np.pi*N_c*r_c*(S_w-1.0)/(F_k+F_d)

#    print 'q_v : ', np.max(q_v), np.min(q_v)
#    print 'q_c : ', np.max(q_c), np.min(q_c)
#    print 'Delta: ',np.max(dt*ConC), np.min(dt*ConC)

    return ConC

    
def rainwater_condensation(rho=None,p=None,T=None,q_r=None,q_c=None,
                           q_v=None,q_sat=None,rho_g=None,q_g=None,dt=4.0,**kwargs):
    
    from scipy.special import gamma

    P_sat=611.2*np.exp(17.62*(T-273.0)/(T-30.0))
    S_w=np.minimum(q_v/q_sat,1.0)

    ilambda_r=(np.maximum((rho_g*q_r)/(np.pi*q_g*rho_w*N_0r),0))**0.25
    F_k=(L_v/(R_v*T)-1.0)*(L_v/(K_a*T))
    F_d=(R_v*T)/(P_sat*D_v)


    print 'rain condensation in: %f %f'%(q_g.min(),rho_g.min())
    ConR=(2.0*np.pi*q_g/rho_g*N_0r*(S_w-1.0)/(F_k+F_d)*
            ilambda_r**2.0
            +0.22*gamma(2.75)*np.sqrt(a_r/nu)
                  *(rho_0/rho_g)**0.25*ilambda_r**(11.0/12.0))
    print 'rain condensation out'

    return ConR

def water_droplet_autoconversion(q_c=None,q_r=None,rho_g=None,dt=4.0,**kwargs):

    q_g=1.0-q_c-q_r

    AutC=k_c*(q_c-q_g/rho_g*a_c)
    return np.maximum(AutC,0.0)


def warm_accreation(q_g=None,q_c=None,q_r=None,rho_g=None,rho_0=1.0,**kwargs):

    from scipy.special import gamma

    print 'accreation in'
    ilambda_r=(np.maximum((rho_g*q_r)/(np.pi*q_g*rho_w*N_0r),0))**0.25
    Acc=np.pi/4.0*N_0r*a_r*np.sqrt(rho_0/rho_g)*gamma(3.5)*ilambda_r**3.5*q_c
    print 'accreation out'

    return Acc


class hot_moist_microphysics_model(microphysics_model):

    prescribed_fields={'T':"InsituTemperature",
                    'rho':"Density",
                    'p':"EOSPressure"}
    forced_fields={'q_v':"WaterVapourFraction",
                   'q_c':"CloudWaterFraction",
                   'q_r':"RainWaterFraction",
                   'e_i':"InternalEnergy"}
    python_fields=[
        ('q_g',get_q_g),
        ('rho_g',get_gas_density),
        ('q_sat',get_q_sat),
        ('q_w',get_q_w),
        ]


    # First build your microphysical source terms:


    ConN=microphysics_forcing(supersaturation_fix,(('q_c',1.0),
                                                   ('q_v',-1.0),
                                                   ('e_i',L_v)),
                              ('q_c','q_v'))
    ConC=microphysics_forcing(cloudwater_condensation,(('q_c',1.0),
                                                       ('q_v',-1.0),
                                                       ('e_i',L_v)),
                               ('q_c','q_v'))
    ConR=microphysics_forcing(rainwater_condensation,(('q_r',1.0),
                                                      ('q_v',-1.0),
                                                      ('e_i',L_v)),
                              ('q_r','q_v'))
    Aut=microphysics_forcing(water_droplet_autoconversion,(('q_r',1.0),
                                                           ('q_c',-1.0)),
                             ('q_r','q_c'))
    Acc=microphysics_forcing(warm_accreation,(('q_r',1.0),
                                              ('q_c',-1.0)),
                             ('q_r','q_c'))


#    forcings=tuple()
#    forcings=(ConC,ConN)
#    forcings=(ConC,ConR,Aut,Acc)
    forcings=(ConN,ConC,ConR,Aut,Acc)

    def __init__(self,state,dt):
        microphysics_model.__init__(self,state, 
                                  self.prescribed_fields,
                                  self.forced_fields,
                                  self.python_fields,
                                  self.forcings,
                                  tuple(),
                                  dt
                                  )

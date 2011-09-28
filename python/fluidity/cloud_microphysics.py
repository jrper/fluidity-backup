#default_physical_parameters={'m':0.018}
import numpy as np
import fluidity_tools

def hot_moist_microphysics(state,parameters):
    dynamic_fields={'T':"InSituTemperature",
                    'rho':"Density",
                    'p':"EOSPressure"}
    active_tracer_fields={'q_v':"WaterVapourFraction",
                          'q_c':"CloudWaterFraction",
                          'q_r':"RainWaterFraction"}


    for field,name in dynamic_fields.items():
        exec('%s=state.scalar_fields["%"][:]'%(field,name))
        exec('dynamic_fields["%s"]=%s'%(name,name))
    for field,name in active_tracer_fields.items():
        exec('%s=state.scalar_fields["%"][:]'%(field,name))
        exec('delta%s=state.scalar_fields["%sSource"]'%(field,name))
        exec('active_tracer_fields["%s"]=%s'%(name,name))

    args=dynamic_fields
    args.update(active_tracer_fields)

    args['q_sat']=get_q_sat(**args)
    ConR=rainwater_condensation(**args)
    ConN=supersaturation_fix(**args)
    ConC=cloud_water_condensation(**args)
    Aut=water_droplet_autocorrelation(**args)
    Acc=accreation(**args)

    nodes=range(q_v.node_count)

    dq_v[nodes]=-ConR-ConC-ConN
    dq_c[nodes]=ConN+ConC+Aut-Acc
    dq_r[nodes]=ConR+Acc
    
        


def get_q_sat(rho=None,p=None,T=None,R_v=None,**kwargs):
    return (p*m)/(rho*R_v*T)

def supersaturation_fix(q_v=None,q_sat=None,dt=None,**kwargs):
    return (q_v-q_sat)/dt


def cloud_water_evaporation(rho=None,p=None,T=None,q_g=None,
                            q_v=None,q_c=None,q_sat=None,**kwargs):
    S_w=q_v/q_sat
    rho_g=rho*q_g
    r_c=10.0e-6
    N_c=(3.0*rho*q_c)/(4.0*rho_w*np.pi*r_c**3)
    F_k=(L_v/(R_v*T-1.0))*(L_v/(K_a*T))
    F_d=(R_v*T)/(P_sat*D_v)

    return 4.0*np.pi*q_g/rho_g*N_c*r_c*(S_w-1.0)/(F_k+F_d)
    
def rainwater_condensation(rho=None,p=None,T=None,q_r=None,
                           q_v=None,q_sat=None,**kwargs):
    
    lambda_r=((np.pi*q_g*rho_w*N_0r)/(rho_g*q_r))**0.25
    rho_g=rho*q_g
    F_k=(L_v/(R_v*T-1.0))*(L_v/(K_a*T))
    F_d=(R_v*T)/(P_sat*D_v)

    return (2.0*np.pi*q_g/rho_g*N_0r*(S_w-1.0)/(F_k+F_d)*
            1.0/lambda_r**2
            +0.22*np.gamma*=(2.75)*np.sqrt(a_r)/eta
                  *(rho_0/rho_g)**0.25*lambda_r**(-11.0/12.0))

def water_droplet_autocorrelation(**kwargs):
    pass


def warm_accreation(**kwargs):
    pass


from math import exp, sqrt, sin, cos, pi
import numpy

p0=1.0e5
cp=1000.0
cv=714.285714
R=cp-cv
theta0=300.0
g=9.81

cp_v=1840.0
R_v=461.5
cv_v= cp_v-R_v
c_w=4100.0
rho_w=1000.0

a=0.75
b=1.0-a

def cv_g(X=None,qv=None,qc=None):
	if X != None:
		return (cv*q_g(X)+(cv_v-cv)*q_v(X))/(q_g(X))
	else:
		return (cv+(cv_v-cv)*qv-cv*qc)/(1.0-qc)

def cv_T(X=None,qv=None,qc=None):
	if X != None:
		return (cv+(cv_v-cv)*q_v(X)
			  +(c_w-cv)*q_c(X)
			  +(c_w-cv)*q_r(X))
	else:
		return cv+(cv_v-cv)*qv

def cp_T(X=None,qv=None,qc=None):
	if X != None:
		return (cp+(cp_v-cp)*q_v(X)
			  +(c_w-cp)*q_c(X)
			  +(c_w-cp)*q_r(X))
	else:
		return cp+(cp_v-cp)*qv

def cp_g(X=None,qv=None,qc=None):
	if X != None:
		return (cp*q_g(X)+(cp_v-cp)*q_v(X))/q_g(X)
	else:
		return (cp+(cp_v-cp)*qv-cp*qc)/(1.0-qc)

def R_g(X=None,qv=None,qc=None):
	return cp_g(X,qv,qc)-cv_g(X,qv,qc)

def ifix(X,rho_g):
	return rho_g*(q_g(X))/(q_g(X)+rho_g*(q_c(X)+q_r(X))/rho_w)


def pert(X,sig=1.0,dx=0,dy=0):
	x_c=500.0
	y_c=350.0
	r=sqrt((X[0] - x_c - dx)**2 + (X[1] - y_c - dy)**2)
	r_c=200.0
	if r < r_c:
		out=sig*(1.0+cos(pi*r/r_c))/2.0
	else:
		out=0.0

	return out

def theta(X):
	theta_bar=theta0
	theta_dash=b*pert(X,3.0,125.0,250.0)

	return theta_bar+theta_dash

def expres(X):
	return 1.0-g*X[1]/(theta0*cp) 


def TT(X):
	return (theta(X)*cp_T(X)/
		((p0/p(X))**(R_g(X)/cp_g(X))*cp_g(X)*q_g(X)
		+c_w*q_r(X)+c_w*q_c(X)))

def T0(X):
	return theta0*((p(X)/p0)**(R/cp))

def rho(X):
	return ifix(X,cv_g(X)*p(X)/(R_g(X)*ie(X)))
#	return cv_g(X)*p(X)/(R_g(X)*ie(X))
#	return cv*p(X)/(R*ie(X))
	
def rho0(X):
	return cv*p(X)/(R*ie0(X))

def p(X):
	return p0*expres(X)**(cp/R)

def ie(X):
	return cv_g(X)*TT(X)#theta(X)*expres(X)

def ie0(X):
	return cv*T0(X)#theta0*expres(X)


def q_v(X):
	return a*pert(X,0.01)

def q_c(X):
	return a*pert(X,0.001)

def q_r(X):
	return a*pert(X,0.0)

def q_g(X):
	return 1.0-q_c(X)-q_r(X)



	




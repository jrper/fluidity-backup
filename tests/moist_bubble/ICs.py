from math import exp, sqrt, sin, cos, pi
import numpy


L_v=2500000.0
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

a=1.0
b=1.0-a

def cv_g(X=None,qv=None,qc=None,qr=0.0):
	if X != None:
		return (cv*q_a(X)+cv_v*q_v(X))/(q_g(X))
	else:
		return (cv+(cv_v-cv)*qv-cv*(qc+qr))/(1.0-qc-qr)

def cv_T(X=None,qv=None,qc=None):
	if X != None:
		return (cv*q_a(X)+cv_v*q_v(X)
			  +c_w*q_c(X)
			  +c_w*q_r(X))
	else:
		return cv+(cv_v-cv)*qv+(c_w-cv)*qc

def cp_T(X=None,qv=None,qc=None):
	if X != None:
		return (cp*q_a(X)+cp_v*q_v(X)
			  +c_w*q_c(X)
			  +c_w*q_r(X))
	else:
		return cp+(cp_v-cp)*qv+(c_w-cp)*qc

def cp_g(X=None,qv=None,qc=None,qr=0.0):
	if X != None:
		return (cp*q_a(X)+cp_v*q_v(X))/q_g(X)
	else:
		return (cp+(cp_v-cp)*qv-cp*(qc+qr))/(1.0-qc-qr)

def R_g(X=None,qv=None,qc=None):
	return cp_g(X,qv,qc)-cv_g(X,qv,qc)

def ifix(X,rho_g):
	return rho_g/(q_g(X)+rho_g*(q_c(X)+q_r(X))/rho_w)

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
#	return (theta(X)*cp_T(X)/
#		((p0/p(X))**(R_g(X)/cp_g(X))*cp_g(X)*q_g(X)
#		+c_w*q_r(X)+c_w*q_c(X)))

	from scipy.optimize import fsolve

	def TtoTheta(T):
#		return (T*((p(X)/p0)**(-R_g(X)/(cp_T(X))))
#			*exp(L_v*(q_c(X)+q_r(X))/(cp_T(X)*T)))-theta(X)


		return T*((p0/p(X))**(R_g(X)/cp_g(X))*cp_g(X)*q_g(X))


	def dTtoTheta(T):
#		return (((p(X)/p0)**(-R_g(X)/(cp_T(X))))
#			*exp(L_v*(q_c(X)+q_r(X))/(cp_T(X)*T))
#			     -L_v*(q_c(X)+q_r(X))/(cp_T(X)*T**2)*
#			     T*((p(X)/p0)**(-R_g(X)/(cp_T(X))))
#			*exp(L_v*(q_c(X)+q_r(X))/(cp_T(X)*T)))

		return ((p0/p(X))**(R_g(X)/cp_g(X))*cp_g(X)*q_g(X))

	T0=theta(X)*(p(X)/p0)**(R_g(X)/(cp_T(X)))
		    
			  

#	T=fsolve(TtoTheta,T0,
#		 fprime=dTtoTheta)

	T= theta(X)*(p(X)/p0)**(R_g(X)/cp_g(X))
	return T

def T0(X):
	return theta0*((p(X)/p0)**(R/cp))

def rho(X):
#	return ifix(X,p(X)/(R_g(X)*TT(X)))
#	return cv_g(X)*p(X)/(R_g(X)*ie(X))
#	return cv*p(X)/(R*ie(X))

	return p(X)/(R_g(X)*q_g(X)*TT(X)+p(X)*(q_c(X)+q_r(X))/rho_w)
	
def rho0(X):
	return p(X)/(R*T0(X))

def p(X):
	return p0*expres(X)**(cp/R)

def ie(X):
	return cv_T(X)*TT(X)#theta(X)*expres(X)

def ie0(X):
	return cv*T0(X)#theta0*expres(X)


def q_v(X):
	return a*pert(X,0.018)

def q_c(X):
	return 0.0+a*pert(X,0.001)

def q_r(X):
	return a*pert(X,0.0)

def q_g(X):
	return 1.0-q_c(X)-q_r(X)

def q_a(X):
	return 1.0-q_c(X)-q_r(X)-q_v(X)



	




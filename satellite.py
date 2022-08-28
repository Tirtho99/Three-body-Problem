import numpy as np
import matplotlib.pyplot as p
def s_1(r,phi,t):
	return np.sqrt(r**2+re**2-2*r*re*np.cos(phi-w*t))**3
def s_2(r,phi,t):
		return np.sqrt(r**2+rm**2-2*r*rm*np.cos(phi-w*t-np.pi))**3
def Rvel(Pr):
	return Pr/mu
def phivel(Pphi,r):
	return Pphi/(mu*r**2)
def pr_dot(Pphi,r,phi,t):
	return (Pphi**2/(mu*r**3))-(Gm1*mu*(r-r1*np.cos(phi-w*t)))/((s_1(r,phi,t)))-(Gm2*mu*(r-r2*np.cos(phi-w*t-np.pi)))/((s_2(r,phi,t)))
def Pphi_dot(r,phi,t):
	return -(Gm1*mu*(r*r1*np.sin(phi-w*t)))/((s_1(r,phi,t)))-(Gm2*mu*(r*r2*np.sin(phi-w*t-np.pi)))/((s_2(r,phi,t)))
Rm=5
mu=1#Mass of thing
w=2*np.pi/(24*27.0)
r1=re=384000.0/(1+Rm)
r2=rm=384000.0*Rm/(1+Rm)
Gm1=w**2*r2*14.7e10
Gm2=Gm1/Rm
phx=(np.pi/180)*104
r=r0=np.sqrt(6400.0**2+r1**2-2*6400*r1*np.cos(-phx+np.pi))
phi=phi0=np.arcsin((6400/r0)*np.sin(-phx+np.pi))
t=t0=0
Pr=Pr0=mu*-37000
Pphi=Pphi0=-(np.pi/180)*150#Zero angular velocity
tf=24*15
h=0.01
X=[r0*np.cos(phi0)]
Y=[r0*np.sin(phi0)]
Xe=[re*np.cos(w*t0)]
Ye=[re*np.sin(w*t0)]
Xm=[rm*np.cos(w*t0+np.pi)]
Ym=[rm*np.sin(w*t0+np.pi)]
dist=[np.sqrt((Xm[-1]-X[-1])**2+(Ym[-1]-Y[-1])**2)]
dist1=[np.sqrt((Xe[-1]-X[-1])**2+(Ye[-1]-Y[-1])**2)]
n=int((tf-t0)/h)
for i in range(0,n):
	if dist[-1]<1700:
		break
	if dist1[-1]<6700 and t >5:
		break
	dr1=Rvel(Pr)*h
	dph1=phivel(Pphi,r)*h
	dPr1= pr_dot(Pphi,r,phi,t)*h
	dPphi1=Pphi_dot(r,phi,t)*h
	dr2=Rvel(Pr+0.5*dPr1)*h
	dph2=phivel(Pphi+0.5*dPphi1,r+0.5*dr1)*h
	dPr2= pr_dot(Pphi+0.5*dPphi1,r+0.5*dr1,phi+0.5*dph1,t+0.5*h)*h
	dPphi2=Pphi_dot(r+0.5*dr1,phi+0.5*dph1,t+0.5*h)*h
	dr3=Rvel(Pr+0.5*dPr2)*h
	dph3=phivel(Pphi+0.5*dPphi2,r+0.5*dr2)*h
	dPr3= pr_dot(Pphi+0.5*dPphi2,r+0.5*dr2,phi+0.5*dph2,t+0.5*h)*h
	dPphi3=Pphi_dot(r+0.5*dr2,phi+0.5*dph2,t+0.5*h)*h
	dr4=Rvel(Pr+dPr3)*h
	dph4=phivel(Pphi+dPphi3,r+dr3)*h
	dPr4= pr_dot(Pphi+dPphi3,r+dr3,phi+dph3,t+h)*h
	dPphi4=Pphi_dot(r+dr3,phi+dph3,t+h)*h
	r+=(1/6.0)*(dr1+2*dr2+2*dr3+dr4)
	phi+=(1/6.0)*(dph1+2*dph2+2*dph3+dph4)
	Pr+=(1/6.0)*(dPr1+2*dPr2+2*dPr3+dPr4)
	Pphi+=(1/6.0)*(dPphi1+2*dPphi2+2*dPphi3+dPphi4)
	if abs(X[-1]-r*np.cos(phi))>1000 or abs(Y[-1]-r*np.sin(phi))>1000:
		n=n-int((t-tf)/h)
		h=h/2.0
		n=n+int((t-tf)/h)
		continue 
	t+=h
	X.append(r*np.cos(phi))
	Y.append(r*np.sin(phi))
	Xe.append(re*np.cos(w*t))
	Ye.append(re*np.sin(w*t))
	Xm.append(rm*np.cos(np.pi+w*t))
	Ym.append(rm*np.sin(np.pi+w*t))
	dist.append(np.sqrt((Xm[-1]-X[-1])**2+(Ym[-1]-Y[-1])**2))
	dist1.append(np.sqrt((Xe[-1]-X[-1])**2+(Ye[-1]-Y[-1])**2))
m=dist.index(np.amin(dist))
p.plot(X[m],Y[m], marker="o",markersize=20, markeredgecolor="red")
p.plot(Xm[m],Ym[m], marker="o",markersize=20, markeredgecolor="red")
p.plot(X,Y)
p.plot(Xm,Ym)
p.plot(Xe,Ye)
p.grid()
p.title("Hit Moob")
#p.axis([-30000,50000,-2000,18000])
p.plot(Xe[-1],Ye[-1], marker="o",markersize=20, markeredgecolor="blue")
p.annotate("Earth",(Xe[-1],Ye[-1]))
p.plot(Xm[-1],Ym[-1], marker="o",markersize=20, markeredgecolor="red")
p.annotate("Moon",(Xm[-1],Ym[-1]))
p.savefig("map1.jpeg")
p.show()


	
	
	
	




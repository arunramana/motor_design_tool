from ...trig import *

def tau_mag(p,wd0,i_d,iq,current_angle):
    #p -> poles pairs
    i_pk = (i_d**2+iq**2)**0.5
    w = wd0
    
    return (3/2)*p*w*i_pk*cos(current_angle)


def tau_rel(p,ld,lq,i_d,iq,current_angle):
    #p -> poles pairs

    i_pk = (i_d**2+iq**2)**0.5

    return (3/4)*p*(i_pk**2)*sin(2*current_angle)*(lq-ld)

def getvd(p,i_d,iq,omega,lq):
   
    rs =20*1e-3
    
    #vd = rs*i_d - p*(omega/9.55)*lq*(1e-6)*iq
    vd = -p*(omega/9.55)*lq*(1e-6)*iq
    return vd

def getvq(p,i_d,iq,omega,ld):
    
    rs = 20*1e-3
    psi = 0.02074340818430611
    
    #vq = rs*iq + p*(omega/9.55)*(ld*1e-6)*i_d+p*(omega/9.55)*psi
    vq = p*(omega/9.55)*(ld*1e-6)*i_d+p*(omega/9.55)*psi
    return vq

def getvoltage(iq,p,i_d,omega,lq,ld,voltage):
    
    return (getvq(p,i_d,iq,omega,ld)**2 + getvd(p,i_d,iq,omega,lq)**2)**0.5 - voltage

def getVolt(omega,p,i_d,iq,lq,ld):
    
    return (getvq(p,i_d,iq,omega,ld)**2 + getvd(p,i_d,iq,omega,lq)**2)**0.5
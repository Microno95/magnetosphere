import numpy as np

def BS_Farris1994(th,D_OB,M,gamma=5./3):
    r = ((gamma-1)*M**2+2)/((gamma+1)*(M**2-1))
    
    if(type(th)==type(np.array([]))):
        RC = radius_of_curvature(th,D_OB)
        return D_OB[1:-1]+RC*0.8*r
    else:
        return D_OB+1.1*r
                           
def BS_GasDyn(D_OB,M,gamma=5./3):
    ''' Gas Dynamic Bow Shock (From Maksimovic/Spreiter) '''
    r = ((gamma-1)*M**2+2)/((gamma+1)*(M**2-1))
    return D_OB*(1+1.1*r)
    

def radius_of_curvature(th,r):
    from numpy import gradient,sqrt,abs
    
    dth = abs(th[1]-th[0])
    
    dr_dth = (r[2:]-r[:-2])/(2*dth)
    d2r_dth2 = (r[2:]-2*r[1:-1]+r[:-2])/dth**2
    
    A = sqrt(r[1:-1]**2 + dr_dth**2)**3
    B = abs( r[1:-1]**2 +2*dr_dth**2 - r[1:-1]*d2r_dth2 )
    
    return A/B

def BS_Jerab05(phi, th, Pd, Ma, B, gamma=5./3):
    ''' The Jerab et. al 2005 bow shock model
    phi: Azimuth from +y toward +z
    th: colatitude from +x
    N: SW number density in cm-3
    v: SW speed in km/s
    B: IMF strength in nT
    '''

    a11 = 0.45
    a22 = 1.
    a33 = 0.8
    a12 = 0.18
    a14 = 46.6
    a24 = -2.2
    a34 = -0.6
    a44 = -618.
    
    def Rav(phi, th):
        x = np.cos(th)
        y = np.sin(th)*np.sin(phi)
        z = np.sin(th)*np.cos(phi)
        
        a = a11*x**2+a22*y**2+a33*z**2 + a12*x*y
        b = a14*x + a24*y + a34*z
        c = a44
        return (-b + np.sqrt(b**2-4*a*c))/(2*a)

    D = 0.937*(0.846+0.042*B)
    C = 91.55
    
    R0 = Rav(0,0)

    A = C*(1.67e-6/Pd)**(1/6)
    m = 1+D*((gamma-1)*Ma**2+2)/((gamma+1)*(Ma**2-1))

    return Rav(phi, th)/R0*A*m
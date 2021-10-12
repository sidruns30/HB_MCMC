import numpy as np
import matplotlib.pyplot as plt

PI = 3.14159265358979323846
G = 6.6743e-8
C = 2.998e10
AU = 1.496e13
MSUN = 1.9885e33
RSUN = 6.955e10
SEC_DAY = 86400.0

''' Some useful functions
'''
pow = lambda x,j: np.power(x,j)
sin = lambda x: np.sin(x)
cos = lambda x: np.cos(x)
tan = lambda x: np.tan(x)
exp = lambda x: np.exp(x)
log = lambda x: np.log10(x)
p10 = lambda x: np.float_power(10, x)
sqrt = lambda x: np.sqrt(x)
asin = lambda x: np.arcsin(x)
atan = lambda x: np.arctan(x)
SQR = lambda x: x*x
CUBE = lambda x: x*x*x
QUAD = lambda x: x*x*x*x

'''Functions to find the overlap (1/3)
'''
def swap (x1, x2):
    temp = x1
    x1 = x2
    x2 = temp
    return x1, x2

''' Functions to find the overlap (2/3)
'''
def A_rh(R, h):
    return R*R*asin(h/R) - h*sqrt(SQR(R) - SQR(h))

''' Functions to find overlap (3/3)
'''
def overlap(r1, r2, d):
    area = 0
    if (r2 > r1): r1, r2 = swap(r1, r2)
    d = abs(d)
    if (d >= r1+r2): area = 0
    if (d < r1-r2): area = PI*r2*r2
    dc = sqrt(r1*r1 - r2*r2)
    h_sq = (SQR(2*d*r1) - SQR(d - r2*r2 + r1*r1)) / SQR(2*d)
    h = sqrt(h_sq)
    if ((d > dc)&(d < (r1+r2))): area = A_rh(r1,h)+A_rh(r2,h)
    if ((d <= dc)&(d >= (r1-r2))): area = PI*r2*r2-(A_rh(r2,h)-A_rh(r1,h))
    return area


''' Trajectory function useful in calculating the lightcurve
'''
def traj(t, x, pos, zdot, rE, theta, rr, ff):
    t = t*SEC_DAY
    # Extract the parameters
    M1 = p10(x[0])
    M2 = p10(x[1])
    P = p10(x[2])*SEC_DAY
    e = x[3]
    inc = x[4]*PI/180
    Omega = x[5]*PI/180
    omega0 = x[6]*PI/180
    T0 = x[7]*SEC_DAY

    Mtot = M1 + M2
    a = pow(G*Mtot*SQR(P)/SQR(2*PI), 1/3)
    M = 2*PI*(t-T0)/P
    M = M % (2*PI)
    EE = M
    sin_M = sin(M)

    if sin_M != 0: EE = M + 0.85*e*sin_M/abs(sin_M)

    # Recursively solve the Kepler's equation
    for i in range(6): EE -= (EE - e*sin(EE) - M)/(1-e*cos(EE))

    r = a*(1-e*cos(EE))
    f = 2*atan(sqrt((1+e)/(1-e))*tan(EE/2))

    cos_Omega = cos(Omega) 
    cos_omega0_f = cos(omega0+f) 
    sin_Omega = sin(Omega) 
    sin_omega0_f = sin(omega0+f) 
    cos_inc = cos(inc) 
    sin_inc = sin(inc) 
    cos_omega0 = cos(omega0)

    XX = r*(cos_Omega*cos_omega0_f - sin_Omega*sin_omega0_f*cos_inc)
    YY = r*(sin_Omega*cos_omega0_f + cos_Omega*sin_omega0_f*cos_inc)
    ZZ = r*sin_omega0_f*sin_inc

    zdot = 1./(sqrt(1-e*e))*(cos_omega0_f + e*cos_omega0)
    pos[0] = XX*(M2/Mtot)
    pos[1] = YY*(M2/Mtot)
    pos[2] = ZZ*(M2/Mtot)
    pos[3] = -XX*(M1/Mtot)
    pos[4] = -YY*(M1/Mtot)
    pos[5] = -ZZ*(M1/Mtot)

    rE = sqrt(4*G*M1/(C*C)*abs(ZZ))
    theta = omega0+f-PI/2
    rr = r
    ff = f
    return zdot, pos, rE, theta, rr, ff

''' Generate and plot the lightcurves taken from likelihood2.c
'''
def lightcurve_model(times, Nt, x):
    # Stored lightcurve
    template = []

    # Helpful constants
    R1 = 1.0 #units of RSUN
    R2 = 1.0 #units of RSUN
    Teff1 = 1.0 #units of Tsun
    Teff2 = 1.0 #units of Tsun
    Flux1 = 1.0 #units of Fsun
    Flux2 = 1.0 #units of Fsun
    Flux_TESS = 1.0 #TESS flux units
    alphabeam = 1.
    alphaev = 1.
    zdot=0.
    rE=0.
    theta=0.
    rr=0.
    ff=0.
    aR=0.
    Am = 0.
    Amag1, Amag2 = np.zeros(Nt), np.zeros(Nt)
    pos = [1,0,0,-1,0,0]
    Tcoeff = [3.74677,0.557556,0.184408,-0.0640800,-0.0359547]
    Rcoeff = [0.00158766,0.921233,-0.155659,-0.0739842,0.0581150]
    itime=0

    # Find a power law fit for temperatures and radii
    logM1 = x[0]
    logM2 = x[1] 
    R1 = 0.0
    Teff1 = 0.0
    R2 = 0.0
    Teff2 = 0.0

    for j in range(4):
        R1 += Rcoeff[j]*pow(logM1,j)
        Teff1 += Tcoeff[j]*pow(logM1,j)
        R2 += Rcoeff[j]*pow(logM2,j)
        Teff2 += Tcoeff[j]*pow(logM2,j)

    # Extract the parameters
    M1 = p10(x[0])
    M2 = p10(x[1])
    P = p10(x[2])*SEC_DAY
    e = x[3]
    inc = x[4]*PI/180
    Omega = x[5]*PI/180
    omega0 = x[6]*PI/180
    T0 = x[7]*SEC_DAY
    Flux_TESS = p10(x[8])
    rr1 = p10(x[9])
    rr2 = p10(x[10])
    print("Period is %f " %P)
    # More parameter scaling
    Mtot = (M1+M2)*MSUN
    a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.)
    #P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot))
    Pdays = P/SEC_DAY
    #printf("%12.5e %12.5e %12.5e\n",a,P,Pdays)
    R1 = pow(10.,R1)*rr1
    R2 = pow(10.,R2)*rr2
    Teff1 = pow(10.,Teff1)/5580.
    Teff2 = pow(10.,Teff2)/5580.

    Flux1 = PI*R1*R1*QUAD(Teff1)
    Flux2 = PI*R2*R2*QUAD(Teff2)

    #Siddhant: taking out all the redundant function calls from the loop
    #convert back to Agnieszka units
    Mtot = M1+M2
    aR = a/RSUN
    Amag_limb=1.

    sin_inc = sin(inc)
    iPdays_to_one_third = pow(Pdays, -1./3)
    iMtot_to_two_thirds = pow(Mtot, -2./3)

    for itime in range(0, Nt, 1):
        t = times[itime]
        zdot, pos, rE, theta, rr, ff = traj(t, x, pos, zdot, rE=0, theta=0, rr=0, ff=0)

        rr = rr/RSUN   #units of RSUN
        cos_ff = cos(ff) 
        cos_2theta = cos(2.*theta) 
         #projected separation in units of RSUN
        d = sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+(pos[4]-pos[1])*(pos[4]-pos[1]))/RSUN 
         #*****************COMPUTE COEFFICIENTS************************/
         # For Amag2 */
        Adoppler = 2.8e-3 * alphabeam * sin_inc * iPdays_to_one_third * iMtot_to_two_thirds * M1 * zdot 
        Aellipse_phi = -alphaev * (M1/M2) * SQR(sin_inc) * cos_2theta * CUBE(R2) / CUBE(rr) 
        Am = (1./9.) * alphaev * (2.+5.*M1/M2) * (2.-3.* SQR(sin_inc)) *CUBE(R2) / CUBE(aR) 
        Aellipse_mean = Am * (e*cos_ff * (3. + 3.*e*cos_ff + e*e*cos_ff*cos_ff)
                                    + 3.*e - 3.*e*e + e*e*e) 
        Aellipse_mean = Aellipse_mean/(CUBE(1.0-e*e)) 
        Amag2[itime] = Flux2 * (Amag_limb + Adoppler + Aellipse_phi + Aellipse_mean) 
        
         #Siddhant: why are we computing the same thing again? I might remove this*/
         # For Amag 1*/
        Adoppler = -1 * Adoppler * M2 / M1   #-2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M2*zdot 
        Aellipse_phi = -alphaev * (M2/M1) * SQR(sin_inc) * cos_2theta * CUBE(R1) / CUBE(rr) 
        Am = (1./9.) * alphaev * (2.+5.*M2/M1) * (2.-3.* SQR(sin_inc)) * CUBE(R1) / CUBE(aR) 
        Aellipse_mean = Am * (e*cos_ff * (3. + 3.*e*cos_ff + e*e*cos_ff*cos_ff)
                                    + 3.*e - 3.*e*e + e*e*e) 
        Aellipse_mean = Aellipse_mean/CUBE(1.0-e*e) 
        Amag1[itime] = Flux1 * (Amag_limb+Adoppler + Aellipse_phi + Aellipse_mean)

        area = overlap(R1,R2,d)
        if (pos[5] > pos[2]): Amag2[itime]-=area*QUAD(Teff2)
        elif (pos[5] < pos[2]): Amag1[itime]-=area*QUAD(Teff1)
        print("FLuxes are: %f \t %f \n" %(Amag1[itime], Amag2[itime]))
        template.append((Amag1[itime]+Amag2[itime])*Flux_TESS)
    return np.array(template)

# Set parameter values to generate the lightcurve
M1 = 1#-1.5, 2
M2 = 1#-1.5, 2
P = -1#-2, 3
e = 0#0, 1
ic = 0#-90, 90
Omega = 0#-180, 180
omega0 = 0#-180, 180
T0 =  0#-1000, 1000
Flux_TESS = 0 #-5, 5
rr1 =  1#-2, 2
rr2 =  1#-2, 2

x = [M1, M2, P, e, ic, Omega, omega0, T0, Flux_TESS, rr1, rr2]

Nt = 100
times = np.linspace(0, 0.1, Nt)
lc = lightcurve_model(times, Nt, x)


plt.figure(figsize=(15,5))
plt.plot(times, lc)
plt.ylim((-5, 5))
plt.savefig("temp_lc.png")
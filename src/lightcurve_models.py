'''
Python script that produces lightcurves based on choice of input parameters and lightcurve models.

The parameters and their units are:
M1:             Mass of star 1 (log Msun)
M2:             Mass of star 2 (log Msun)
P:              Period (log days)
e:              Eccentricity
inc:            Inclination (degrees)
Omega:          Longitude of ascending node (radians)
omega0:         Argument of periapse (radians)
T0:             Phase of the lightcurve (days)
rr1:            Radius scaling factor for star 1; R1 ~ R1_i x 10^(rr1)
rr2:            Radius scaling factor for star 2

There are two lightcurve models based on "Morris et. al 1985" and "Engel et. al 2020"
Morris et. al model includes doppler and ellipsoidal variations whereas the Engel model
includes reflection as well.
Eclipse correction to the flux are added after each model is applied
A common trajectory function is used to calculate the trajectory

Task list:
    Compute the trajectory
    Calculate Teff_i and R_i from M_i and find the effective flux for each star
    Apply the lightcurve model
    Add in the eclipse correction
'''

# Set parameters here
parameters = {
    'M1'        : -9.24270e-01,
    'M2'        : -7.47407e-01,    
    'P'         : 7.96050e-01 ,
    'e'         : 3.17503e-01,
    'inc'       : 1.49535e+00*180/3.141,
    'Omega'     :  2.23680e+00 ,    
    'omega0'    : -2.48707e+00 ,
    'T0'        : 6.56090e+01 ,
    'rr1'       : 6.44932e-01 ,
    'rr2'       : 8.65771e-01 ,
}

# Choice of model "engel" or "morris"
model = "morris"

# Plot individual contributions (ellipsoidal variation, eclipses, etc.) in the model
plot_comp = True
# Plot overall lightcurves from both the models
plot_both = True
# Name to save the lightcurve plot as
lc_name = model + "_lc.png"

################################## MAIN CODE STARTS BELOW #######################################
#################################################################################################
import numpy as np
from numpy import sin, cos, tan, log10, power, exp, fmod, fabs,sqrt
from numpy import arctan as atan
from numpy import arcsin as asin
import matplotlib.pyplot as plt
import os
from matplotlib import gridspec


# Problem constants
PI = 3.14159265358979323846
G = 6.6743e-8
C = 2.998e10
AU = 1.496e13
MSUN = 1.9885e33
RSUN = 6.955e10
SEC_DAY = 86400.0
NPARS = 11

# Useful functions
SQR = lambda x: x*x
CUBE = lambda x: x*x*x
QUAD = lambda x: x*x*x*x

'''
Temperature and radius function:
Parameters: M (mass of star Msun), rr (radius scaling factor)
Returns:    Teff (K), R (radius)
'''
def get_TR(M, rr):
    logM = log10(M)
    Tcoeff = [3.74677,0.557556,0.184408,-0.0640800,-0.0359547]
    Rcoeff = [0.00158766,0.921233,-0.155659,-0.0739842,0.0581150]
    R = 0.
    Teff = 0.
    for j in range(5):
        # R is in log Rsun
        R += Rcoeff[j]*pow(logM,j)
        Teff += Tcoeff[j]*pow(logM,j)
    # [Units are Rsun and K respectively]
    R = pow(10.,R)*rr
    Teff = pow(10.,Teff)/5580.
    return Teff, R 


'''
Trajectory function: computes positions of the two stars in cgs
Parameters: current time (days), model parameters
Returns: Cartersian and polar positions of the two stars and the radial velocity
'''
def traj(t, pars):
    t = t*SEC_DAY
    M1 = pow(10.,pars[0])*MSUN
    M2 = pow(10.,pars[1])*MSUN
    P = pow(10.,pars[2])*SEC_DAY
    e = pars[3]
    inc = pars[4]*(PI/180.)
    Omega = pars[5]
    omega0 = pars[6]
    T0 = pars[7]*SEC_DAY    
    Mtot = M1+M2
    a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.)
    M = 2.*PI * (t-T0)/P 
    M = fmod(M,2*PI)
    EE = M
    sin_M = np.sin(M)

    if(sin_M != 0.0):   EE = M + 0.85*e*sin_M/fabs(sin_M)

    for i in range(4):  EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE))

    r = a*(1-e*cos(EE))
    f = 2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.))

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

    X1 = XX*(M2/Mtot)
    Y1 = YY*(M2/Mtot)
    Z1 = ZZ*(M2/Mtot)
    X2 = -XX*(M1/Mtot)
    Y2 = -YY*(M1/Mtot)
    Z2 = -ZZ*(M1/Mtot)
    rr = r
    ff = f

    zdot = 1 / sqrt(1 - SQR(e)) * (cos_omega0_f + e*cos_omega0)
    vrad_over_M = sqrt(G/(Mtot * a)) * sin(inc) * zdot

    return X1, Y1, Z1, X2, Y2, Z2, rr, ff, zdot, vrad_over_M

'''
Ecclipse area function: computers the overlap between the two stars
Parameters: R1, R2 (radii of stars), X1, X2, Y1, Y2 (cartesian postions of the stars)
Returns: Overlap area 
'''
def eclipse_area(R1, R2, X1, X2, Y1, Y2):
    # Overlap function borrowed from likelihood2.c
    d = sqrt(SQR(X2-X1) + SQR(Y2-Y1))/RSUN
    
    if (R2 > R1):
        temp_ = R1
        R1 = R2
        R2 = temp_
    
    area = 0.
    d = fabs(d)
    dc = sqrt(R1*R1-R2*R2)
    # Now find the observed overlapping area between the two stars
    if (d >= (R1+R2)): area = 0.
    if (d < (R1-R2)): area = PI*R2*R2

    if ((d > dc)&(d < (R1+R2))):
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d)
        h = sqrt(h_sq)

        Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h)
        Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h)
        area = Arh1 + Arh2
        

    if ((d <= dc)&(d >= (R1-R2))):
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d)
        h = sqrt(h_sq)
        Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h)
        Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h)
        area = PI*R2*R2-(Arh1 + Arh2)
        
    return area

'''
Engel model
Parameters: model parameter array and true anomaly array
Returns: Total and partial fluxes of the two stars
'''
class engel:
    def __init__(self, pars, nu_arr, vrad_over_m_arr):
        self.pars = pars
        # Extract the paramters
        self.logM1 = pars[0]
        self.logM2 = pars[1]
        # Period in seconds
        self.P = pow(10., pars[2])*SEC_DAY
        self.Pdays = pow(10., pars[2])
        self.e = pars[3]
        self.inc = pars[4]*(PI/180)
        self.Omega = pars[5]
        self.omega0 = pars[6]
        self.T0 = pars[7]*SEC_DAY
        self.rr1 = pow(10., pars[8])
        self.rr2 = pow(10., pars[9])
        self.M1 = pow(10., self.logM1)
        self.M2 = pow(10., self.logM2)

        self.nu_arr = nu_arr
        self.vrad_over_m = vrad_over_m_arr

    '''
    Modified based on the results from Loeb, Gaudi 2003, ApJ
    A = alpha_beam 4 vrad / c 
    Engel et. al may have a typo!
    '''
    def beaming(self, M2, vrad_over_M, alpha_beam, nu, M1, Pdays, inc, omega0, e):
        # The commented code is an equivalent definition
        #M2 *= MSUN
        #return -alpha_beam * M2 * vrad_over_M * 4 / C
        q = M2 / M1
        return -2830 * alpha_beam * q / pow(1+q, 2/3) * pow(M1, 1/3) * pow(Pdays, -1/3) * sin(inc) * cos(omega0 + nu) / sqrt(1 - SQR(e)) * pow(10, -6)


    def ellipsoidal(self, mu, tau, P, M1, M2, e, inc, omega0, nu, R1, a, **kwargs):
        '''fixed bugs in: alpha_11'''
        alpha_11 = 15 * mu * (2 + tau) / (32 * (3 - mu))
        alpha_21 = 3 * (15 + mu) * (1 + tau) / (20 * (3 - mu))
        alpha_2b1 = 15 * (1 - mu) * (3 + tau) / (64 * (3 - mu))
        alpha_01 = alpha_21 / 9
        alpha_0b1 = 3 * alpha_2b1 / 20
        alpha_31 = 5 * alpha_11 / 3
        alpha_41 = 7*alpha_2b1 / 4

        beta = (1 + e * cos(nu)) / (1 - SQR(e))
        q = M2 / M1
        # Angular velocity at the periapse
        Prot = P * pow(1 - e, 3./2)

        ppm = 1.e-6
        # Flag to include const terms in flux
        enable_const = 1
        tot_flux = 0.
        '''fixed bugs in: CONS2'''
        if (enable_const == 1):
            # Note that onlt the first term is constant
            CONS1 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * CUBE(R1) / (M1 * SQR(Prot)) * ppm
            CONS2 = 13435 * 3 * alpha_01 * (2 - 3*SQR(sin(inc))) * q / (1 + q) * CUBE(beta*R1) / (M1 * SQR(P)) * ppm
            CONS3 = (759 * alpha_0b1 * (8 - 40*SQR(sin(inc)) + 35*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3) 
                    * pow(P, -10./3) * pow(beta * R1, 5)) * ppm

            tot_flux += (CONS1 + CONS2 + CONS3) 
        

        S1 = (3194 * alpha_11 * (4 * sin(inc) - 5 * CUBE(sin(inc))) * pow(M1, -4./3) * q / pow(1+q, 4./3) *  pow(P, -8./3)
                * QUAD(beta * R1) * sin(omega0 + nu)) * ppm
        C2_1 = 13435 * alpha_21 * SQR(sin(inc)) * q / (1 + q) * CUBE(beta * R1) / (M1 * SQR(P)) * cos(2*(omega0 + nu)) * ppm
        C2_2 =  (759 * alpha_2b1 * (6*SQR(sin(inc)) - 7*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3)
                * (beta * R1) * QUAD(beta * R1) * cos(2*(omega0 + nu))) * ppm
        S3 = (3194 * alpha_31 * CUBE(sin(inc)) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3) * QUAD(beta * R1) 
                * sin(3*(omega0 + nu))) * ppm
        C4 = (759 * alpha_41 * QUAD(sin(inc)) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3) * (beta * R1) *
                QUAD(beta * R1) * cos(4*(omega0 + nu))) * ppm

        tot_flux += (S1 + C2_1 + C2_2 + S3 + C4)
        
        if ("components" in kwargs and enable_const == 1): return tot_flux, CONS1, CONS2, CONS3, S1, C2_1, C2_2, S3, C4
        else: return tot_flux

    def reflection(self, P, M1, M2, e, inc, omega0, nu, R2, alpha_ref1):
        q = M2 / M1
        beta = (1 + e * cos(nu)) / (1 - SQR(e))
        ppm = 1.e-6

        fac1 = pow(1 + q, -2./3)
        fac2 = pow(M1, -2./3)
        fac3 = pow(P, -4./3)
        fac4 = SQR(beta *R2)
        fac5 = 0.64 - sin(inc) * sin(omega0 + nu) + 0.18 * SQR(sin(inc)) * (1 - cos(2*(omega0 + nu)))

        return 56514 * alpha_ref1 * fac1 * fac2 * fac3 * fac4 * fac5 * ppm

    def generate_lightcurve(self, ellip_args=False):
        # First calculate the effective temperatures and radii
        Teff1, R1 = get_TR(self.M1, self.rr1)
        Teff2, R2 = get_TR(self.M2, self.rr2)

        Flux1 = PI*SQR(R1)*QUAD(Teff1)
        Flux2 = PI*SQR(R2)*QUAD(Teff2)
        Ftot = Flux1 + Flux2
        Flux1 /= Ftot
        Flux2 /= Ftot

        Nt = len(self.nu_arr)
        Amag1, Amag2, Abeam, Aellip, Aref = map(np.ones, [Nt, Nt, Nt, Nt, Nt])

        Mtot = (self.M1+self.M2)*MSUN
        a = pow(G*Mtot*self.P*self.P/(4.0*PI*PI),1./3.)
        ar = a / RSUN

        alpha_beam = 1.     #Beaming
        alpha_ref = 0.1      #Reflection 
        mu = .16           #Limb darkening
        tau = .344            #Gravity darkening
    
        CONS1_arr, CONS2_arr, CONS3_arr, S1_arr, C2_1_arr, C2_2_arr, S3_arr, C4_arr = map(np.zeros, [Nt, Nt, Nt, Nt, Nt, Nt, Nt, Nt])

        for i, nu in enumerate(self.nu_arr):
            beam1 = self.beaming(self.M2, self.vrad_over_m[i], alpha_beam, nu, self.M1, self.Pdays, self.inc, self.omega0, self.e)
            ellip1, CONS1, CONS2, CONS3, S1, C2_1, C2_2, S3, C4 = self.ellipsoidal(.835, .435, self.Pdays, self.M1, self.M2, self.e, self.inc, self.omega0, nu, R1, ar, components=True)
            ref1 = self.reflection(self.Pdays, self.M1, self.M2, self.e, self.inc, self.omega0, nu, R2, .48)

            # Add in the individua ellipsoidal constributions
            CONS1_arr[i] += Flux1 *CONS1
            CONS2_arr[i] += Flux1 *CONS2
            CONS2_arr[i] += Flux1 *CONS3
            S1_arr[i] += Flux1 *S1
            C2_1_arr[i] += Flux1 *C2_1
            C2_2_arr[i] += Flux1 *C2_2
            S3_arr[i] += Flux1 *S3
            C4_arr[i] += Flux1 *C4

            beam2 = self.beaming(self.M1, self.vrad_over_m[i], alpha_beam, nu, self.M2, self.Pdays, self.inc, self.omega0+PI, self.e)
            ellip2, CONS1, CONS2, CONS3, S1, C2_1, C2_2, S3, C4  = self.ellipsoidal(.94, .56, self.Pdays, self.M2, self.M1, self.e, self.inc, (self.omega0+PI), nu, R2, ar, components=True)
            ref2 = self.reflection(self.Pdays, self.M2, self.M1, self.e, self.inc, (self.omega0+PI), nu, R1, .03)

            # Add in the individua ellipsoidal constributions
            CONS1_arr[i] += Flux2 *CONS1
            CONS2_arr[i] += Flux2 *CONS2
            CONS2_arr[i] += Flux2 *CONS3
            S1_arr[i] += Flux2 *S1
            C2_1_arr[i] += Flux2 *C2_1
            C2_2_arr[i] += Flux2 *C2_2
            S3_arr[i] += Flux2 *S3
            C4_arr[i] += Flux2 *C4

            Amag1[i] = Flux1 * (1 + beam1 + ellip1 + ref1)
            Amag2[i] = Flux2 * (1 + beam2 + ellip2 + ref2)
            Abeam[i] = Flux1 * beam1 + Flux2 * beam2
            Aellip[i] = Flux1 * ellip1 + Flux2 * ellip2
            Aref[i] = Flux1 * ref1 + Flux2 * ref2

        if ellip_args:
            return Amag1, Amag2, Abeam, Aellip, Aref, CONS1_arr, CONS2_arr, CONS3_arr, S1_arr, C2_1_arr, C2_2_arr, S3_arr, C4_arr
        else:
            return Amag1, Amag2, Abeam, Aellip, Aref

'''
Morris model
Parameters: model parameter array, radial separation array, zdot array and true anomaly array
Returns: Total and partial fluxes of the two stars
'''
class morris:
    def __init__(self, pars, rr_arr, zdot_arr, ff_arr):
        self.pars = pars
        # Extract the paramters
        self.logM1 = pars[0]
        self.logM2 = pars[1]
        # Period in seconds
        self.P = pow(10., pars[2])*SEC_DAY
        self.Pdays = pow(10., pars[2])
        self.e = pars[3]
        self.inc = pars[4]*(PI/180)
        self.Omega = pars[5]
        self.omega0 = pars[6]
        self.T0 = pars[7]*SEC_DAY
        self.rr1 = pow(10., pars[8])
        self.rr2 = pow(10., pars[9])
        self.M1 = pow(10., self.logM1)
        self.M2 = pow(10., self.logM2)
        
        self.rr_arr = rr_arr
        self.zdot_arr = zdot_arr
        self.ff_arr = ff_arr

    def doppler(self, M1, M2, inc, Pdays, zdot):
        alpha_beam = 1
        Mtot = M1 + M2
        return -2.8e-3*alpha_beam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M2*zdot

    def ellipse_phi(self, M1, M2, inc, theta, R1, rr):
        alpha_ev = 1
        # Radial separation in solar radii
        rr /= RSUN
        return -alpha_ev*(M2/M1)*sin(inc)*sin(inc)*cos(2.*theta)*R1*R1*R1/(rr*rr*rr)

    '''
    Mean distortion on the 
    '''
    def ellipse_mean(self, M1, M2, inc, R1, aR, e, ff):
        alpha_ev = 1
        # How much star gets distorted
        Am = (1./9.)*alpha_ev*(2.+5.*M2/M1)*(2.-3.*sin(inc)*sin(inc))*R1*R1*R1/(aR*aR*aR)
        # How mean star rotates in the orbit
        Am *= (e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e)

        Am /= (1.0-e*e)*(1.0-e*e)*(1.0-e*e)
        return  Am

    def generate_lightcurve(self):
        # First calculate the effective temperatures and radii
        Teff1, R1 = get_TR(self.M1, self.rr1)
        Teff2, R2 = get_TR(self.M2, self.rr2)

        Flux1 = PI*SQR(R1)*QUAD(Teff1)
        Flux2 = PI*SQR(R2)*QUAD(Teff2)
        Ftot = Flux1 + Flux2
        Flux1 /= Ftot
        Flux2 /= Ftot

        Nt = len(self.rr_arr)
        Amag1, Amag2, Adop, Aellip_var, Aellip_mean = map(np.ones, [Nt, Nt, Nt, Nt, Nt])

        Mtot = (self.M1+self.M2)*MSUN
        a = pow(G*Mtot*self.P*self.P/(4.0*PI*PI),1./3.)
        ar = a / RSUN

        for i in range(Nt):
            theta = self.ff_arr[i] + self.omega0 - PI/2
            dop1 = self.doppler(self.M1, self.M2, self.inc, self.Pdays, self.zdot_arr[i])
            ellip_phi1 = self.ellipse_phi(self.M1, self.M2, self.inc, theta, R1, self.rr_arr[i])
            ellip_mean1 = self.ellipse_mean(self.M1, self.M2, self.inc, R1, ar, self.e, self.ff_arr[i])

            dop2 = -self.doppler(self.M2, self.M1, self.inc, self.Pdays, self.zdot_arr[i])
            ellip_phi2 = self.ellipse_phi(self.M2, self.M1, self.inc, theta, R2, self.rr_arr[i])
            ellip_mean2 = self.ellipse_mean(self.M2, self.M1, self.inc, R2, ar, self.e, self.ff_arr[i])

            Amag1[i] = Flux1 * (1 + dop1 + ellip_phi1 + ellip_mean1)
            Amag2[i] = Flux2 * (1 + dop2 + ellip_phi2 + ellip_mean2)
            Adop[i] = Flux1 * dop1 + Flux2 * dop2
            Aellip_var[i] = Flux1 * ellip_phi1 + Flux2 * ellip_phi2
            Aellip_mean[i] = Flux1 * ellip_mean1 + Flux2 * ellip_mean2

        return Amag1, Amag2, Adop, Aellip_var, Aellip_mean

'''
Full lightcurve function:
Parameters: times (time array), pars (model parameters), model_type ("engel" or "morris"), 
            plot_comp (choice to return individual comp)
Returns: array(s) of the lightcurve 
'''
def compute_lc(times, parameters, model_type, plot_comp, **kwargs):
    # convert dict of parameters to lists
    pars = [vals for _, vals in parameters.items()]
    Nt = len(times)

    X1arr = np.zeros(Nt)
    Y1arr = np.zeros(Nt)
    Z1arr = np.zeros(Nt)
    X2arr = np.zeros(Nt)
    Y2arr = np.zeros(Nt)
    Z2arr = np.zeros(Nt)
    rr_arr = np.zeros(Nt)
    nu_arr = np.zeros(Nt)
    zdot_arr = np.zeros(Nt)
    vrad_over_m_arr = np.zeros(Nt)
    Aecl = np.zeros(Nt)

    #Compute the trajectory
    for i in range(Nt): 
        (X1arr[i], Y1arr[i], Z1arr[i], X2arr[i], Y2arr[i],
         Z2arr[i], rr_arr[i], nu_arr[i], zdot_arr[i], vrad_over_m_arr[i]) = traj(times[i], pars)

    # Get the model contributions
    if model_type == "engel":
        Amag1, Amag2, Abeam, Aellip, Aref = engel(pars, nu_arr, vrad_over_m_arr).generate_lightcurve()
        
    elif model_type == "morris":
        Amag1, Amag2, Adop, Aellip_var, Aellip_mean = morris(pars, rr_arr, zdot_arr, nu_arr).generate_lightcurve()

    else: raise ValueError("Incorrect model type")

    # Eclipse contribution
    M1 = power(10., pars[0])
    M2 = power(10., pars[1])
    rr1 = power(10., pars[8])
    rr2 = power(10., pars[9])

    Teff1, R1 = get_TR(M1, rr1)
    Teff2, R2 = get_TR(M2, rr2)

    Flux1 = PI*SQR(R1)*QUAD(Teff1)
    Flux2 = PI*SQR(R2)*QUAD(Teff2)
    Ftot = Flux1 + Flux2
    Flux1 /= Ftot
    Flux2 /= Ftot

    for i in range(Nt):
        area = eclipse_area(R1, R2, X1arr[i], X2arr[i], Y1arr[i], Y2arr[i])
        if (Z2arr[i] > Z1arr[i]): 
            #print(area, R1, R2)
            Amag2[i] -= area * Flux2 / (PI*SQR(R2))
            Aecl[i] = -area #/ (PI*SQR(R2)) 
        elif (Z2arr[i] < Z1arr[i]): 
            Amag1[i] -= area * Flux1 / (PI*SQR(R1))
            Aecl[i] = - area  #/ (PI*SQR(R1)) 

    # Full lightcurve
    Afull = Amag1 + Amag2

    #print(Aecl)

    #Plot the trajectory
    if "plot_traj" in kwargs:
        # Get the parameters and collect arrays to plot
        Emag1, Emag2, Ebeam, Eellip, Eref, CONS1_arr, CONS2_arr, CONS3_arr, S1_arr, C2_1_arr, C2_2_arr, S3_arr, C4_arr = engel(pars, nu_arr, vrad_over_m_arr).generate_lightcurve(ellip_args=True)
        Mmag1, Mmag2, Mdop, Mellip_var, Mellip_mean = morris(pars, rr_arr, zdot_arr, nu_arr).generate_lightcurve()
        inc = parameters["inc"]
        ecc = parameters["e"]
        dirname = "../figures/Trajectory_inc=%3.1e_e=%3.1e" %(inc, ecc) 
        #try:os.mkdir(dirname)
        #except: pass
        X1arr /= RSUN
        Y1arr /= RSUN
        Z1arr /= RSUN
        X2arr /= RSUN
        Y2arr /= RSUN
        Z2arr /= RSUN

        # Load numpy file
        #data = np.loadtxt("output_lc.txt")

        for index in range(0, 1000, 1000):
            # Set figures and axes
            fig = plt.figure(figsize=(20,10), constrained_layout=True)
            gs = gridspec.GridSpec(nrows=10, ncols=10)
            ax1 = fig.add_subplot(gs[:, :5], projection="3d")
            plt.axis('off')
            plt.grid(b=None)
            ax1.set_xlabel("$X [R_{\circ}]$")
            ax1.set_ylabel("$Y [R_{\circ}]$")
            ax1.set_title("Orbital Plane")
            # Hide grid lines
            ax1.grid(False)
            # Hide axes ticks
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_zticks([])
            ax1.view_init(elev=90, azim=0)
            ax1.set_xlim((-20,20))
            ax1.set_ylim((-20,20))
            ax1.set_zlim((-20,20))

            ax2 = fig.add_subplot(gs[0:2, 5:])
            ax2.set_ylabel("$\Delta F / F$")
            ax2.set_title("Doppler Variation")

            ax3 = fig.add_subplot(gs[2:4, 5:])
            ax3.set_ylabel("$\Delta F / F$")
            ax3.set_title("Ellipsoidal Mean Variation")

            ax4 = fig.add_subplot(gs[4:6, 5:])
            ax4.set_ylabel("$\Delta F / F$")
            ax4.set_title("Ellipsoidal Phi Variation")

            ax5 = fig.add_subplot(gs[6:8, 5:])
            ax5.set_ylabel("$RSUN$")
            ax5.set_title("Radial and angular separation")
            ax7 = ax5.twinx()
            ax7.set_ylabel("$\\nu$")

            ax6 = fig.add_subplot(gs[8:, 5:])
            ax6.set_ylabel("$\Delta F$")
            ax6.set_xlabel("$T (days)$")
            ax6.set_title("Full lightcurve")

            plt.subplots_adjust(hspace=1.2)
            plt.suptitle("Period: %2.1e [days]; Eccentricity: %2.1e; Inclination %2.1e $^{\circ}$" %(10**parameters["P"], parameters["e"], parameters["inc"]))

            ax1.plot(X1arr, Y1arr, Z1arr)
            ax1.plot(X2arr, Y2arr, Z2arr)
            
            # Make stars
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)


            xs1 = X1arr[index] + R1 * np.outer(np.cos(u), np.sin(v))
            ys1 = Y1arr[index] + R1  * np.outer(np.sin(u), np.sin(v))
            zs1 = Z1arr[index] + R1  * np.outer(np.ones(np.size(u)), np.cos(v))
            # Plot the star surface
  
            xs2 = X2arr[index] + R2  * np.outer(np.cos(u), np.sin(v))
            ys2 = Y2arr[index] + R2  * np.outer(np.sin(u), np.sin(v))
            zs2 = Z2arr[index] + R2  * np.outer(np.ones(np.size(u)), np.cos(v))
            # Plot the star surface
            if (Z2arr[index] > Z1arr[index]):
                ax1.plot_surface(xs1, ys1, zs1, color='b')
                ax1.plot_surface(xs2, ys2, zs2, color='r')
            else:
                ax1.plot_surface(xs2, ys2, zs2, color='r')
                ax1.plot_surface(xs1, ys1, zs1, color='b')

            ax2.plot(times, Ebeam, "o", label="Engel")
            ax2.plot(times, Mdop, "o",label="Morris")
            #ax2.plot(data[:,0], data[:,2], ".", label="C model")
            ax2.plot(times[index], Ebeam[index], "x")
            ax2.plot(times[index], Mdop[index], "x")

            
            ax3.plot(times, (CONS1_arr + CONS2_arr + CONS3_arr), "o", label="Engel")
            ax3.plot(times, Mellip_mean, "o", label="Morris")
            ax3.plot(times[index], (CONS1_arr + CONS2_arr + CONS3_arr)[index], "x")
            ax3.plot(times[index], Mellip_mean[index], "x")
            
            ax4.plot(times, (C2_1_arr + C2_2_arr + S1_arr + S3_arr + C4_arr +CONS1_arr + CONS2_arr + CONS3_arr), "o", label="Engel")
            ax4.plot(times, Mellip_var, "o", label="Morris")
            #ax4.plot(data[:,0], data[:,3], ".", label="C model")
            ax4.plot(times[index], (C2_1_arr + C2_2_arr + S1_arr + S3_arr + C4_arr)[index], "x")
            ax4.plot(times[index], Mellip_var[index], "x")
            
            '''Plotting other things besides reflection'''
            ax5.plot(times, (X1arr-X2arr)/RSUN, 'r', times,  (Y1arr-Y2arr)/RSUN, "b")
            #ax7.plot(times, X2arr/RSUN, "y", times, Y2arr/RSUN, "g")
            #ax7.plot(times, nu_arr, "r-", label="Engel")
            #ax7.plot(times[index], nu_arr[index]/RSUN, "rx")
            #ax5.plot(data[:,0], data[:,4], ".", label="C model")

            Engel_full = Ebeam + Eellip + Aecl + Eref
            Engel_full -= np.median(Engel_full)
            Engel_full += 1

            Morris_full = Mdop + Mellip_var + Mellip_mean + Aecl
            Morris_full -= np.median(Morris_full)
            Morris_full += 1

            #data[:,1] -= np.median(data[:,1])
            #data[:,1] += 1

            #ax6.plot(times, Engel_full, "-", label="Engel")
            #ax6.plot(times, Morris_full, "o", label="Morris")
            #ax6.plot(data[:,0], data[:,1], ".", label="C model")
            ax6.plot(times, Aecl, "-")
            ax6.plot(times[index], Aecl[index], "-")
            #ax6.set_ylim((.8, 1.01))
            #ax6.set_xlim((1.2, 1.8))

            for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:    ax.legend(prop={"size":8})
            plt.savefig("test.png")
            #plt.savefig(dirname + "/_%03d.png" % (index//10))
            plt.close()

    # Normalize everything
    norm_factor = 1 / np.median(Afull)

    Afull *= norm_factor
    Aecl *= norm_factor

    if plot_comp == True:
        if model_type == "engel":
            Abeam *= norm_factor
            Aellip *= norm_factor
            Aref *= norm_factor

            return Afull, Abeam, Aellip, Aref, Aecl
        
        elif model_type == "morris":
            Adop *= norm_factor
            Aellip_var*= norm_factor
            Aellip_mean *= norm_factor
            
            return Afull, Adop, Aellip_var, Aellip_mean, Aecl

    else: return Afull




#########################################################################################

#Generate a time array for 3 periods
period = power(10., parameters["P"])
times = np.linspace(0, 3*period, 1000)
compute_lc(times, parameters, model_type="engel", plot_comp=True, plot_traj=True)

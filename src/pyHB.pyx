# cython: language_level = 3

cimport cython
#cdef void calc_light_curve(double* times, double Nt, double*pars, double *template);
import numpy as np

cimport likelihood2
cimport likelihood3
import sys
import traceback

def lightcurve2(times,inpars):
  logM1, logM2, logP_day, e, inc, omega, omega0, T0_day,log_rad1_rescale,log_rad2_rescale,logTanom,blend_frac,logFluxTESS=inpars
  cdef double pars[12]
  #need to convert angles from rad to deg here
  radeg=180/np.pi
  pars[:]=[  logM1, logM2, logP_day, e, inc*radeg, omega*radeg,omega0*radeg, T0_day,logFluxTESS,log_rad1_rescale,log_rad2_rescale, logTanom ]
  times=np.array(times)
  cdef int Nt=len(times)
  cdef double[:] ctimes = times
  cdef double[:] ctemplate=np.empty(Nt, dtype=np.double)
  #compute_lightcurve(&ctimes[0],&cAmags[0],Nt,&pars[0]);
  likelihood2.calc_light_curve(&ctimes[0],Nt,&pars[0],&ctemplate[0]);
  template=np.array(ctemplate)
  #return template+10**(log_blendFlux+logFluxTESS)
  Flux_TESS = 10.**logFluxTESS
  return 1*blend_frac + template*(1-blend_frac)


def lightcurve3(times,inpars):
  logM1, logM2, logP_day, e, inc, omega0, T0_day, log_rad1_rescale, log_rad2_rescale, mu1, tau1, mu2, tau2, alprefl1, alprefl2, ln_beam_resc_1, ln_beam_resc_2, ln_alp_Teff_1, ln_alp_Teff_2, blend_frac, flux_tune=inpars
  cdef double pars[20]
  #need to convert angles from rad to deg here
  radeg=180/np.pi
  pars[:]=[  logM1, logM2, logP_day, e, inc, 0, omega0, T0_day, log_rad1_rescale,log_rad2_rescale, mu1, tau1, mu2, tau2, alprefl1, alprefl2, np.exp(ln_beam_resc_1), np.exp(ln_beam_resc_2), np.exp(ln_alp_Teff_1), np.exp(ln_alp_Teff_2) ]
  times=np.array(times)
  cdef int Nt=len(times)
  cdef double[:] ctimes = times
  cdef double[:] ctemplate=np.empty(Nt, dtype=np.double)
  '''
    c-code param definitions:

    // Extract the paramters
    double logM1 = pars[0];
    double logM2 = pars[1];
    // Period in seconds
    double P = pow(10., pars[2])*SEC_DAY;
    double Pdays = pow(10., pars[2]);
    double e = pars[3];
    double inc = pars[4];
    double Omega = pars[5];
    double omega0 = pars[6];
    double T0 = pars[7]*SEC_DAY;
    double rr1 = pow(10., pars[8]);
    double rr2 = pow(10., pars[9]);
    // Limb and gravity darkening coefficients respectively
    mu_1 = pars[10];
    tau_1 = pars[11];
    mu_2 = pars[12];
    tau_2 = pars[13];
    // Reflection coefficients
    alpha_ref_1 = pars[14];
    alpha_ref_2 = pars[15];
  '''
  likelihood3.calc_light_curve(&ctimes[0],Nt,&pars[0],&ctemplate[0]);
  template=np.array(ctemplate)
  return ( 1*blend_frac + template*(1-blend_frac) ) * flux_tune

class parspace:
  def __init__(self, *args):
    self.names=[]
    self.mins=[]
    self.maxs=[]
    if len(args)%2 !=0:raise ValueError("parspace:Constructor requires arguments in pattern ('name1',[min,max],'name2',[min,max],...)")
    for i in range(int(len(args)/2)):
      self.names.append(args[i*2])
      self.mins.append(args[i*2+1][0])
      self.maxs.append(args[i*2+1][1])
    self.maxs=np.array(self.maxs)
    self.mins=np.array(self.mins)
    self.N=len(self.names)
    self.live=np.array([True]*self.N)
    self.pinvals=[None]*self.N
    self.idx={name:i for i,name in enumerate(self.names)}
    self.Nlive=self.N
  def reset_range(self,name,minmax):
    i=self.idx[name]
    if not self.live[i]:
      if self.pinvals[i]<minmax[0] or self.pinvals[i]>minmax[1]:
        raise ValueError('pinned value is not within range')
    self.mins[i]=minmax[0]
    self.maxs[i]=minmax[1]
  def pin(self,name,value):
    i=self.idx[name]
    if(value<self.mins[i] or value>self.maxs[i]):
      print('parspace.pin: Value '+name+' = '+str(value)+'  out of range ['+str(self.mins[i])+','+str(self.maxs[i])+']')
      return False
    if(self.live[i]):self.Nlive-=1
    self.live[i]=False
    self.pinvals[i]=value
    return True
  def get_pars(self,livevals):
    parvals=np.array(self.pinvals)
    parvals[self.live]=livevals
    return parvals
  def live_ranges(self):
    return np.vstack([self.mins[self.live],self.maxs[self.live]]).T
  def live_names(self):
    return [name for i,name in enumerate(self.names) if self.live[i]]
  def draw_live(self):
    #print('Nlive',self.Nlive)
    pars=np.random.rand(self.Nlive)
    #print('s1',pars.shape,self.mins.shape,self.maxs.shape,self.live.shape)
    mn=self.mins[self.live]
    mx=self.maxs[self.live]
    #print('s2',mn.shape,mx.shape)
    return pars*(mx-mn)+mn
  def out_of_bounds(self,pars):
    return not np.all(np.logical_and(np.array(pars)>=self.mins, np.array(pars)<=self.maxs))
    

sp2=parspace(
  'logM1', [ -1.5, 2.0 ],
  'logM2', [ -1.5, 2.0 ],
  'logP', [ -2.0, 3.0 ],
  'e', [ 0, 1 ],
  'inc', [ 0, np.pi ],
  'Omega', [ -np.pi, np.pi ],
  'Omega0', [ -np.pi, np.pi ],
  'T0', [ -1000, 1000 ],
  'log_rad1_resc', [ -2, 2 ],
  'log_rad2_resc', [ -2, 2 ],
  #'log_rad1_resc', [ -0.25, 1.25 ],
  #'log_rad2_resc', [ -0.25, 1.25 ],
  'logTanom', [ -0.5, 0.5 ],
  'blend_frac', [ 0, 1.0 ],
  'logFluxTESS', [ -10.0, 10.0 ],
  'ln_noise_resc', [ -0.2, 0.2 ]
)

sp3=parspace(
  'logM1', [ -1.5, 2.0 ],
  'logM2', [ -1.5, 2.0 ],
  'logP', [ -2.0, 3.0 ],
  'e', [ 0, 1 ],
  'inc', [ 0, np.pi ],
  'omega0', [ -np.pi, np.pi ],
  'T0', [ -1000, 1000 ],
  'log_rad1_resc', [ -2, 2 ],
  'log_rad2_resc', [ -2, 2 ],
  'mu_1', [ 0.12, 0.20 ],
  'tau_1', [ 0.30, 0.38 ],
  'mu_2', [ 0.12, 0.20 ],
  'tau_2', [ 0.30, 0.38 ],
  'alpha_ref_1', [0.8,1.2],
  'alpha_ref_2', [0.8,1.2],
  'ln_beam_resc_1', [-0.1,0.1],
  'ln_beam_resc_2', [-0.1,0.1],
  'ln_alp_Teff_1', [-0.1,0.1],
  'ln_alp_Teff_2', [-0.1,0.1],
  'blend_frac', [ 0.0, 1.0 ],
  'flux_tune', [ 0.99, 1.01 ],
  'ln_noise_resc', [ -0.2, 0.2 ]
)

def likelihood(times,fluxes,errs,pars,lctype=3):
  minlike=-1e18
  ln_noise_resc=pars[-1]
  pars=pars[:-1]
  try:
    if lctype==2:
      modelfluxes=lightcurve2(times, pars)
    elif lctype==3:
      modelfluxes=lightcurve3(times, pars)
    else:
      raise ValueError('Unknown light curve model type.')
    noise_resc=np.exp(ln_noise_resc)
    sigmas=errs*noise_resc
    llike = - np.sum(((fluxes-modelfluxes)/sigmas)**2)/2 - len(errs)*ln_noise_resc # - sum(np.log(2*np.pi*errs**2))/2 #last term is/would be constant.
  except:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    print('likelihood exception:')
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
    llike=minlike
  if not llike>minlike:llike=minlike
  #print(llike)
  return llike

def test_roche_lobe(pars,verbose=False):
    logM1=pars[0]
    logM2=pars[1]
    M1=10**logM1
    M2=10**logM2
    q=M2/M1
    P=pars[2]
    e=pars[3]
    logRadResc1=pars[9]
    logRadResc2=pars[10]
    Rcoeff = [0.00158766,0.921233,-0.155659,-0.0739842,0.0581150]
    R1=R2=0
    for j in range(len(Rcoeff)):
      R1 += Rcoeff[j]*pow(logM1,j)
      R2 += Rcoeff[j]*pow(logM2,j)
    R1=10**(R1+logRadResc1)
    R2=10**(R2+logRadResc2)
    if q<=1:
      Rsec=R2
      Rpri=R1
    else:
      Rsec=R1
      Rpri=R2
    a=4.208278*((M1+M2)*P**2)**(1/3)
    #We artificially symmetrize the Hill radius to somewhat fairly
    #extend the approximate formula applicable to small mass ratio q
    #to a reasonable and symmetric result near equal masses.
    #We get the small mass ratio correctly when q is away from 1
    #with a maximum of 1/2 the minumum separation when the q=1
    Hillfac=((q+2/3+1/q)*3)**(-1/3)
    rperi=a*(1-e)
    aHillsec=rperi*Hillfac
    aHillpri=rperi*(1-Hillfac)#This is a hack version of the test for the primary.
    if verbose:print('Roche lobe test: Rsec, RHillsec, Rpri, RHillpri, :',Rsec,aHillsec,Rpri,aHillpri)
    test=max([Rsec/aHillsec,Rpri/aHillpri])
    return test


import numpy as np
import matplotlib.pyplot as plt
import growth_tools as gt

BoxSize = 5e-3 # simulation box size in Mpc or Mpc/h
GridSize = 512 # GridSize^3 particles in simulation
zinit = 1e6    # starting redshift
zfinal = 100   # ending redshift

gt.use_h_units = True # Mpc/h instead of Mpc -- True for Unigrid, Delta2particles, and (unless changed) Gadget
gt.neutrinos = True # adjust radiation-dominated growth to account for neutrino free-streaming -- True for Delta2particles
gt.alt_fourier_convention = True # divide power spectrum by 8pi^3 -- True for Unigrid input
use_h_units_for_file = False # False for Unigrid input (???)

# k values for power spectrum table
kmin = 1e0
kmax = 1e7
knum = 5000

# power spectrum source
use_file    = True  # use input file (calculated with CAMB Sources) instead of analytic calculation (which neglects baryons)
file_z      = 1000  # redshift of input file
inputFileName = 'planck15_1000_ex.txt' # path to input file

# path into which to save modified power spectrum (as power.dat) and snapshot redshifts (as outputs.txt)
path = '.'

# power spectrum modification
modtype = 'spike'
def modify(k,Pk): # function that takes arrays k, P(k) and modifies P(k)
  if modtype == 'bend':
    kb = 1.
    power = .5
    Pk[k>kb] *= (k[k>kb]/k[k>kb][0])**power
  elif modtype == 'step':
    kb = 1e2
    mult = 64
    Pk[k>kb] *= mult
  elif modtype == 'spike':
    spikek = 30000.
    spikeheight = 625.
    spikelnw = .05
    Pk *= (1.+spikeheight/spikelnw*np.exp(-.5*(np.log(k/spikek)/spikelnw)**2))
  elif modtype == 'plat':
    km = 9000.
    height = 330.
    width = 10.
    Pk[(k>km/width**.5)&(k<km*width**.5)] *= height

# plotting options
plot_Pk     = False # plot power spectrum
plot_Pkdim  = True # plot dimensionless power spectrum
plot_sigma  = False # plot sigmaR
plot_dndlnM = False # plot Press-Schechter halo mass function
plot_evo    = False # plot linear growth
kplot       = 2*np.pi/BoxSize  # k-mode to plot
plot_hor    = False # plot horizon
compare_rel = False # show relativistic power spectra (these neglect neutrinos and baryons)
plot_prim   = False # plot primordial curvature power spectrum
plot_xi     = False # plot correlation function (which is roughly the typical peak shape)
log_xi      = False # plot correlation function in log space (otherwise linear)
plot_snaps  = False # plot snapshot spacing function
analyze_power = False # show statistics of power enhancement
compare_file = False # compare final power spectrum to file:
compareFileName = 'ExtendedICs/planck15_500_ex.dat'

# options for snapshot timing file
snapshots_file = True # generate outputs.txt specifying snapshot times
n_snapshots = 100 # number of snapshots
def outputsFunc(a): # snapshots evenly spaced in this function
  return np.log(gt.D(a,2*np.pi/BoxSize))

# --- end user input

def calcs():
  global MinHaloMass, kminBox, kmaxReq, MassPerParticle, MassInBox
  # calculate parameters
  MassInBox = BoxSize**3*gt.rhoM # Msun pr Msun/h^3
  if gt.use_h_units:
    MassInBox /= gt.HubbleParam**2 # Msun/h
  MassPerParticle = MassInBox/(GridSize)**3
  MinHaloMass = 20*MassPerParticle
  kminBox = 2.*np.pi/BoxSize
  kmaxReq = np.sqrt(3.)*np.pi*GridSize/BoxSize
calcs()

# generate power spectrum
if __name__ == '__main__':
  if plot_Pk:
    plt.figure(1)
    plt.clf()
def calc_pk():
  global k, lnk, Pk, PkFin, PkBase
  k = np.logspace(np.log10(kmin),np.log10(kmax),knum)
  #PkFin = Plate(k,a(zfinal))*power_spectrum_mult
  #Pk = 2*np.pi**2/k**3*Dearly(a(zinit),k)**2
  if use_file:
    from scipy.interpolate import InterpolatedUnivariateSpline as interp
    matter = 3 # 1=CDM, 2=B, 3=M
    data = np.loadtxt(inputFileName)
    PkFile = data[:,matter]
    kFile = data[:,0]
    PkFile *= gt.HubbleParam**3
    kFile /= gt.HubbleParam
    PkFile *= 8*np.pi**3
    PkFile_interp = interp(kFile,PkFile)
    PkFin = PkFile_interp(k)*gt.Dnewt(gt.a(zfinal),k)**2/gt.Dnewt(gt.a(file_z),k)**2
    Pk = PkFile_interp(k)*gt.Dnewt(gt.a(zinit),k)**2/gt.Dnewt(gt.a(file_z),k)**2
  else:
    PkFin = 2*np.pi**2/k**3*gt.Dnewt(gt.a(zfinal),k)**2
    Pk = 2*np.pi**2/k**3*gt.Dnewt(gt.a(zinit),k)**2
  if compare_rel:
    PkNG = 2*np.pi**2/k**3*gt.Dearly(gt.a(zinit),k)**2
    PkSG = 2*np.pi**2/k**3*gt.deltachi(gt.a(zinit),k)**2
  if gt.alt_fourier_convention:
    Pk /= (8*np.pi**3)
    PkFin /= (8*np.pi**3)
    if compare_rel:
      PkNG /= (8*np.pi**3)
      PkSG /= (8*np.pi**3)
  if compare_file:
    matter = 3 # 1=CDM, 2=B, 3=M
    data = np.loadtxt(compareFileName)
    PkFile = data[:,matter]
    kFile = data[:,0]
    if gt.use_h_units and not use_h_units_for_file:
      PkFile *= gt.HubbleParam**3
      kFile /= gt.HubbleParam
    elif not gt.use_h_units and use_h_units_for_file:
      PkFile /= gt.HubbleParam**3
      kFile *= gt.HubbleParam
    if not gt.alt_fourier_convention:
      PkFile *= 8*np.pi**3
    if __name__ == '__main__':
      if plot_Pk:
        plt.figure(1)
        plt.loglog(kFile,PkFile,'g')
  #inds = (Pk>0)
  #Pk = Pk[inds]
  #PkFin = PkFin[inds]
  #k = k[inds]
  lnk = np.log(k)
  # plot
  if __name__ == '__main__':
    if plot_Pk:
      plt.figure(1)
      plt.loglog(k,PkFin,'r--')
      plt.loglog(k,Pk,'b--')
    if plot_Pkdim:
      plt.figure(2)
      plt.clf()
      if compare_rel:
        plt.loglog(k,gt.Pdim(k,PkNG),'g--',lw=3,alpha=.5)
        plt.loglog(k,gt.Pdim(k,PkSG),'y--',lw=3,alpha=.7)
      plt.loglog(k,gt.Pdim(k,PkFin),'r--')
      plt.loglog(k,gt.Pdim(k,Pk),'b--')
      if compare_file:
        plt.figure(2)
        plt.loglog(kFile,gt.Pdim(kFile,PkFile),'g')
  
  # modify power spectrum
  PkBase = Pk.copy()
  modify(k,Pk)
  modify(k,PkFin)
  if compare_rel:
    modify(k,PkNG)
    modify(k,PkSG)
  if __name__ == '__main__':
    # plot
    if plot_Pk:
      plt.figure(1)
      plt.loglog(k,PkFin,'r',label="final (z=%.0f)"%(zfinal))
      plt.loglog(k,Pk,'b',label="initial (z=%.0f)"%(zinit))
    if plot_Pkdim:
      plt.figure(2)
      if compare_rel:
        plt.loglog(k,gt.Pdim(k,PkNG),'g',lw=3,alpha=.5,label="newtonian gauge")
        plt.loglog(k,gt.Pdim(k,PkSG),'y',lw=3,alpha=.7,label="synchronous gauge")
      plt.loglog(k,gt.Pdim(k,PkFin),'r',label="final (z=%.0f)"%(zfinal))
      plt.loglog(k,gt.Pdim(k,Pk),'b',label="initial (z=%.0f)"%(zinit))
    # save modification
    PkSave = Pk.copy()
    PkFinSave = PkFin.copy()
    kSave = k.copy()
    if gt.use_h_units and not use_h_units_for_file:
      PkSave /= gt.HubbleParam**3
      PkFinSave /= gt.HubbleParam**3
      kSave *= gt.HubbleParam
    elif not gt.use_h_units and use_h_units_for_file:
      PkSave *= gt.HubbleParam**3
      PkFinSave *= gt.HubbleParam**3
      kSave /= gt.HubbleParam
    np.savetxt(path+'power.dat',np.stack(np.stack((kSave,PkSave,PkSave,PkSave),axis=1)),fmt='%.15e',delimiter=' ')
    np.savetxt(path+'power_final.dat',np.stack(np.stack((kSave,PkFinSave,PkFinSave,PkFinSave),axis=1)),fmt='%.15e',delimiter=' ')

    # plot labels
    if plot_Pk:
      plt.figure(1)
      plt.loglog(np.full(2,kminBox),np.asarray([Pk[Pk>0].min(),PkFin.max()]),'k--')
      plt.loglog(np.full(2,kmaxReq),np.asarray([Pk[Pk>0].min(),PkFin.max()]),'k--')
      plt.xlabel(r'$k$ (%s)'%gt.unitTexMpc1)
      plt.ylabel(r'$P(k)$ (%s)'%gt.unitTexMpc3)
      plt.legend(loc=3)
    if plot_Pkdim:
      plt.figure(2)
      plt.loglog(np.full(2,kminBox),np.asarray([(gt.Pdim(k[Pk>0],Pk[Pk>0])).min(),(gt.Pdim(k,PkFin)).max()]),'k--')
      plt.loglog(np.full(2,kmaxReq),np.asarray([(gt.Pdim(k[Pk>0],Pk[Pk>0])).min(),(gt.Pdim(k,PkFin)).max()]),'k--')
      plt.ylim(gt.A_s/10,(gt.Pdim(k,PkFin)).max())
      plt.xlabel(r'$k$ (%s)'%gt.unitTexMpc1)
      plt.ylabel(r'$\mathcal{P}(k)$')
      #plt.ylim(1e-12,(k**3*PkFin/(2*np.pi**2)).max())
      plt.legend(loc=2)
calc_pk()

# Press-Schechter
def RfromM(M): #M in Msun
  return ((3.*M)/(4*np.pi*gt.rhoM))**(1./3) #Mpc
  
def flarge(y):
  return 3/y**3*(np.sin(y)-y*np.cos(y))

def fsmall(y):
  return 1 - y**2/10. + y**4/280. - y**6/15120. + y**8/1330560. - y**10/172972800.

def f(y):
  return fsmall(y)*(y<=0.3)+flarge(y)*(y>0.3)
  
def SigmaR(R,z,P,box=False):
  ind = (k <= kmaxReq)&(k>=1./gt.horizon(z))
  if box == True:
    ind = ind&(k>=kminBox)#&(k<=kmaxReq)
  integrand = k[ind]**3*P[ind]*f(k[ind]*R)**2
  integral = np.sum((lnk[ind][1:]-lnk[ind][:-1])*(integrand[:-1]+integrand[1:])/2.) #trapezoidal rule
  if gt.alt_fourier_convention:
    integral *= 8*np.pi**3
  return np.sqrt(1/(2*np.pi**2)*integral)

def SigmaM(M,z,P,box=False):
  ind = (k <= kmaxReq)&(k>=1./gt.horizon(z))
  if box == True:
    ind = ind&(k>=kminBox)
  integrand = k[ind]**3*P[ind]*f(k[ind]*RfromM(M))**2
  integral = np.sum((lnk[ind][1:]-lnk[ind][:-1])*(integrand[:-1]+integrand[1:])/2.) #trapezoidal rule
  if gt.alt_fourier_convention:
    integral *= 8*np.pi**3
  return np.sqrt(1/(2*np.pi**2)*integral)

def fdfdlnylarge(y):
  return flarge(y)*y*(3*np.sin(y)/y**2 - 9*(-y*np.cos(y) + np.sin(y))/y**4)
  
def fdfdlnysmall(y):
  return -y**2/5. + 6.*y**4/175. - 4.*y**6/1575. + 8.*y**8/72765. - y**10/315315.

def fdfdlny(y):
  return fdfdlnysmall(y)*(y<=0.16)+fdfdlnylarge(y)*(y>0.16)

def dlnSigmadlnM(M,P,box=False):
  ind = (k <= kmaxReq)
  if box == True:
    ind = ind&(k>=kminBox)
  integrand = k[ind]**3*P[ind]*fdfdlny(k[ind]*RfromM(M))
  integral = np.sum((lnk[ind][1:]-lnk[ind][:-1])*(integrand[:-1]+integrand[1:])/2.) #trapezoidal rule
  if not gt.alt_fourier_convention:
    integral /= 8*np.pi**3
  return 4*np.pi/(3.*SigmaM(M,zinit,P,box)**2)*integral

def dndlnM(M,z,P,box=False,logM=False):
  if logM:
    M = np.exp(M)
  sigma = SigmaM(M,z,P,box)
  return np.sqrt(2/np.pi)*gt.rhoM*gt.delCrit/(M*sigma)*np.exp(-gt.delCrit**2/(2*sigma**2))*np.abs(dlnSigmadlnM(M,P,box))

def Nhalos(minM,maxM,z,P,box=True):
  from scipy.integrate import quad as integral
  return integral(dndlnM,np.log(minM),np.log(maxM),args=(z,P,box,True))[0]*BoxSize**3

def MdndlnM(M,z,P,box=False,logM=False):
  if logM:
    M = np.exp(M)
  sigma = SigmaM(M,z,P,box)
  return np.sqrt(2/np.pi)*gt.rhoM*gt.delCrit/(sigma)*np.exp(-gt.delCrit**2/(2*sigma**2))*np.abs(dlnSigmadlnM(M,P,box))

def Mhalos(minM,maxM,z,P,box=True):
  from scipy.integrate import quad as integral
  return integral(MdndlnM,np.log(minM),np.log(maxM),args=(z,P,box,True))[0]*BoxSize**3

# formatting
def latex_float(f):
  float_str = "{0:.2g}".format(f)
  if "e" in float_str:
    base, exponent = float_str.split("e")
    return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
  else:
    return float_str

# snapshots file
    
def generateOutputs(scaleFile,zInitial,zFinal,numSnapShots):
  from scipy.optimize import brentq as solve
  aInitial = gt.a(zInitial)
  aFinal = gt.a(zFinal)
  fInitial = outputsFunc(aInitial)
  fFinal = outputsFunc(aFinal)
  fw = open(scaleFile,"w")
  for i in range(numSnapShots):
    if i == 0:
      acurr = aInitial
    elif i == numSnapShots-1:
      acurr = aFinal
    else:
      f = fInitial+(fFinal-fInitial)*i/float(numSnapShots-1)
      def func(a):
        return outputsFunc(a)-f
      acurr = solve(func,aInitial,aFinal)
    fw.write(str(acurr) + "\n")
  fw.close()

if __name__ == '__main__':
  if snapshots_file:
    generateOutputs(path+'outputs.txt',zinit,zfinal,n_snapshots)
  if snapshots_file and plot_snaps:
    snapshots = np.loadtxt(path+'outputs.txt')
    for i in range(snapshots.size):
      print("snapshot %03d: z = %f"%(i,gt.z(snapshots[i])))
    plt.figure(7)
    plt.clf()
    aas = np.logspace(np.log10(gt.a(zinit)),np.log10(gt.a(zfinal)),1000)
    plt.semilogx(aas,outputsFunc(aas))
    plt.plot([gt.aeq,gt.aeq],[outputsFunc(aas).min(),outputsFunc(aas).max()],'c--')
    plt.ylim(outputsFunc(aas).min(),outputsFunc(aas).max())
    plt.xlim(gt.a(zinit),gt.a(zfinal))
    plt.xlabel(r'$a$')
    plt.ylabel(r'snapshot spacing function')
    for i in range(10):
      val = outputsFunc(aas).min() + i*(outputsFunc(aas).max() - outputsFunc(aas).min())/10.
      plt.plot([gt.a(zinit),gt.a(zfinal)],[val,val],'k',alpha=.3)
  # correlation function and average peak
  r = 2*np.pi/k
  #lnr = np.log(r)
  kr = k.reshape((1,-1))*r.reshape((-1,1))
  xi_int = gt.Pdim(k,Pk).reshape((1,-1))*np.sin(kr)/kr
  xi = np.sum((lnk[1:].reshape((1,-1))-lnk[:-1].reshape((1,-1)))*(xi_int[:,:-1]+xi_int[:,1:])/2.,axis=1)
  #dxi_int = gt.Pdim(k,Pk).reshape((1,-1))*k.reshape((1,-1))*(np.cos(kr)/kr-np.sin(kr)/kr**2)
  #dxi = np.sum((lnk[1:].reshape((1,-1))-lnk[:-1].reshape((1,-1)))*(dxi_int[:,:-1]+dxi_int[:,1:])/2.,axis=1)
  #del2xi_int = -k.reshape((1,-1))**2*xi_int
  #del2xi = np.sum((lnk[1:].reshape((1,-1))-lnk[:-1].reshape((1,-1)))*(del2xi_int[:,:-1]+del2xi_int[:,1:])/2.,axis=1)
  if plot_xi:
    plt.figure(8)
    plt.plot(r,xi)
    #plt.plot(r,dxi,'g')
    #plt.plot(r,del2xi,'y')
    plt.xlabel(r'$r$ (%s)'%gt.unitTexMpc)
    plt.ylabel(r'$\xi(r)$')
    plt.xlim((0,.5*BoxSize))
    #plt.yscale('log')
    if log_xi:
      plt.xlim((.5*BoxSize/GridSize,.5*BoxSize))
      plt.xscale('log')
      plt.yscale('log')
    plt.tight_layout()
    xi0 = xi[r>BoxSize/GridSize].max()
    inds = (r>BoxSize/GridSize)&(r<BoxSize*.2)
    xi_avg_int = xi[inds]*4.*np.pi*r[inds]**3
    lnr = np.log(r[inds])
    xi_avg = np.sum((-lnr[1:]+lnr[:-1])*.5*(xi_avg_int[:-1]+xi_avg_int[1:]))/(4./3*np.pi*r[inds].max()**3)
    print("mean xi = %lg"%xi_avg)
      
  # results
  print("Box size = %g %s, initial redshift = %f"%(BoxSize,gt.unitMpc,zinit))
  print("    (1/aH = %g %s, T = %g eV)"%(gt.horizon(zinit),gt.unitMpc,gt.T0/gt.a(zinit)))
  print("    Box is %f of horizon (should be smaller than horizon)"%(BoxSize/gt.horizon(zinit)))
  print("")
  
  print("Check initial and final redshifts:")
  SigmaRInitGrid = SigmaR(0.62*BoxSize/GridSize,zinit,Pk)
  SigmaRBoxInitGrid = SigmaR(0.62*BoxSize/GridSize,zinit,Pk,True)
  SigmaRInitBox = SigmaR(BoxSize,zinit,Pk)
  SigmaRBoxInitBox = SigmaR(BoxSize,zinit,Pk,True)
  print("z = %.0f, R = 0.62 BoxSize/GridSize = %g %s:"%(zinit,0.62*BoxSize/GridSize,gt.unitMpc))
  print("    SigmaR    = %g (should be <= 0.2 and match)"%(SigmaRInitGrid))
  print("    SigmaRBox = %g"%(SigmaRBoxInitGrid))
  print("z = %.0f, R = 0.62 BoxSize/GridSize:"%(zfinal))
  print("    SigmaR    = %g (should be > 1 and match)"%(SigmaR(0.62*BoxSize/GridSize,zfinal,PkFin)))
  print("    SigmaRBox = %g"%(SigmaR(0.62*BoxSize/GridSize,zfinal,PkFin,True)))
  print("z = %.0f, R = BoxSize:"%(zinit))
  print("    SigmaR    = %g (should be << 1)"%(SigmaRInitBox))
  print("    SigmaRBox = %g"%(SigmaRBoxInitBox))
  print("z = %.0f, R = BoxSize:"%(zfinal))
  print("    SigmaR    = %g (should be < 0.1)"%(SigmaR(BoxSize,zfinal,PkFin)))
  print("    SigmaRBox = %g"%(SigmaR(BoxSize,zfinal,PkFin,True)))
  print("")
  print("z = %.0f, M = MassPerParticle:"%(zinit))
  print("    SigmaM    = %g (should be < 0.3 and match)"%(SigmaM(MassPerParticle,zinit,Pk)))
  print("    SigmaMBox = %g"%(SigmaM(MassPerParticle,zinit,Pk,True)))                 
  print("z = %.0f, M = MinHaloMass:"%(zfinal))
  print("    SigmaM    = %g (should be > 1 and match)"%(SigmaM(MinHaloMass,zfinal,PkFin)))
  print("    SigmaMBox = %g"%(SigmaM(MinHaloMass,zfinal,PkFin,True)))
  print("z = %.0f, M = 0.01*MassInBox:"%(zfinal))
  print("    SigmaM    = %g (should be close when > 1)"%(SigmaM(0.01*MassInBox,zfinal,PkFin)))
  print("    SigmaMBox = %g"%(SigmaM(0.01*MassInBox,zfinal,PkFin,True)))
  print("")
  print("z = %.0f, M = MinHaloMass:"%(zfinal))
  print("    dndlnM    = %g (should match)"%(dndlnM(MinHaloMass,zfinal,PkFin)))
  print("    dndlnMBox = %g"%(dndlnM(MinHaloMass,zfinal,PkFin,True)))
  print("")
  print("Press-Schechter analysis at z = %.0f"%zfinal)
  print("    N halos (>=20 particles) = %g"%(Nhalos(MinHaloMass, 0.1*MassInBox, zfinal,PkFin)))
  print("    N large halos (>=200 particles) = %g"%(Nhalos(200*MassPerParticle, 0.1*MassInBox, zfinal,PkFin)))
  print("    N small halos (>=2 particles) = %g"%(Nhalos(2*MassPerParticle, 0.1*MassInBox, zfinal,PkFin)))
  print("    fraction of mass in halos (>=20 particles) = %g"%(Mhalos(MinHaloMass, MassInBox, zfinal,PkFin)/MassInBox))
  #mtotal = MassInBox
  #print("(no box) fraction of mass in halos (>=20 particles) = %g"%(Mhalos(MinHaloMass, mtotal, zfinal,PkFin,False)/mtotal))
  
  '''
  import BBKS_calcs
  PkFin2 = PkFin.copy()
  if gt.alt_fourier_convention:
    PkFin2 *= (2*np.pi)**3
  PkFin2[k>spikek*np.exp(spikelnw)] = 0
  PkFin2[k<spikek/np.exp(spikelnw)] = 0
  bbks = BBKS_calcs.Cosmology(k,PkFin2)
  print("")
  print("N halos (BBKS) = %g"%(BoxSize**3*bbks.n(bbks.nu(gt.delCrit))))
  PR = gt.R0(k)**2
  modify(k,PR)
  PR[k>spikek*np.exp(spikelnw)] = 0
  PR[k<spikek/np.exp(spikelnw)] = 0
  print("spike area %g"%(np.sum(PR)*np.mean(np.diff(lnk))))
  '''
  
  # more analysis
  if analyze_power:
    from scipy.interpolate import interp1d as interp
    ind = (k <= kmaxReq)&(k >= kminBox)
    integrand = k[ind]**3*PkBase[ind]
    integral = np.cumsum((lnk[ind][1:]-lnk[ind][:-1])*(integrand[:-1]+integrand[1:])/2.) #trapezoidal rule
    if gt.alt_fourier_convention:
        integral *= 8*np.pi**3
    powerBase = 1/(2*np.pi**2)*integral
    integrand = k[ind]**3*Pk[ind]
    integral = np.cumsum((lnk[ind][1:]-lnk[ind][:-1])*(integrand[:-1]+integrand[1:])/2.) #trapezoidal rule
    if gt.alt_fourier_convention:
        integral *= 8*np.pi**3
    powerMod = 1/(2*np.pi**2)*integral
    plt.figure(9)
    plt.clf()
    kP = .5*(k[ind][1:]+k[ind][:-1])
    #plt.semilogx(kP,powerMod-powerBase)
    fracPower = (powerMod-powerBase)/(powerMod[-1]-powerBase[-1])
    #plt.figure(10)
    plt.semilogx(kP,fracPower)
    kFracPow = interp(fracPower,kP)
    plt.xlim((kFracPow(.001),kFracPow(.999)))
    plt.xlabel(r'$k$ $(%s)$'%gt.unitMpc1)
    plt.ylabel(r'fraction of excess power in $k^\prime>k$')
    ratioPow = lambda x: kFracPow(.5+.5*x)/kFracPow(.5-.5*x)
    xs = np.linspace(0,.99,100)
    #plt.plot(xs,np.log(ratioPow(xs)))
    #plt.xlim((0,1))
    #plt.xlabel('fraction of excess power')
    #plt.ylabel('number of e-folds in k')
    xprint = [.5,.75,.9,.95,.99]
    print('')
    for p in xprint:
      print('%.2f excess power in %.2f e-folds'%(p,np.log(ratioPow(p))))
  
  if plot_sigma:
    nrad = 100
    radii = np.logspace(np.log10(0.62*BoxSize/GridSize),np.log10(BoxSize),nrad)
    plt.figure(3)
    plt.clf()
    sigmaradi = np.zeros(nrad)
    sigmaradf = np.zeros(nrad)
    for i in range(nrad):
      sigmaradi[i] = SigmaR(radii[i],zinit,Pk,True)
      sigmaradf[i] = SigmaR(radii[i],zfinal,PkFin,True)
    plt.loglog(radii,sigmaradi,'b-',label="z = %.0f, box"%(zinit))
    plt.loglog(radii,sigmaradf,'r-',label="z = %.0f, box"%(zfinal))
    for i in range(nrad):
      sigmaradi[i] = SigmaR(radii[i],zinit,Pk,False)
      sigmaradf[i] = SigmaR(radii[i],zfinal,PkFin,False)
    plt.loglog(radii,sigmaradi,'b--',label="z = %.0f"%(zinit))
    plt.loglog(radii,sigmaradf,'r--',label="z = %.0f"%(zfinal))
    #plt.loglog(np.full(2,0.62*BoxSize/GridSize),np.asarray([sigmaradi.min(),sigmaradf.max()]),'k--')
    #plt.loglog(np.full(2,BoxSize),np.asarray([sigmaradi.min(),sigmaradf.max()]),'k--')
    plt.ylabel(r'$\sigma_R$')
    plt.xlabel(r'$R$ (%s)'%gt.unitTexMpc)
    plt.legend(loc=3)
  
  if plot_dndlnM:
    nlnM = 100
    Mvals = np.logspace(np.log10(MassPerParticle),np.log10(0.1*MassInBox),nlnM)
    dndlnMvals = np.zeros(nlnM)
    for i in range(nlnM):
      dndlnMvals[i] = dndlnM(Mvals[i],zfinal,PkFin,True)
    plt.figure(4)
    plt.clf()
    plt.loglog(Mvals/MassPerParticle,dndlnMvals,label="z = %d"%(zfinal))
    plt.ylabel(r'$\frac{dn}{d\ln M}$')
    plt.xlabel(r'number of particles')
    plt.ylim((.001,np.max(dndlnMvals)))
    plt.legend()
      
  if plot_evo:
    plt.figure(5)
    plt.clf()
    zH = (-gt.OmegaM+np.sqrt(gt.OmegaM**2+4.*gt.OmegaR*(kplot/gt.H0_in_Mpc)**2))/(2.*gt.OmegaR)-1.
    aas = np.logspace(np.log10(gt.a(zH))-.2,np.log10(gt.a(zfinal)),200)
    #aas = np.logspace(np.log10(gt.a(zinit))-.5,np.log10(gt.a(zinit))+.5,200)
    norm = 1.686/gt.D(gt.a(zfinal),kplot)
    delta1s = gt.D(aas,kplot)*norm
    delta4s = gt.Dearly(aas,kplot)*norm
    delta2s = np.zeros(aas.size)
    delta3s = np.zeros(aas.size)
    for i in range(aas.size):
      delta2s[i] = gt.deltachiM(aas[i],kplot)*norm
      delta3s[i] = gt.deltachi(aas[i],kplot)*norm
    plt.plot(aas,delta1s,'b',label='subhorizon soln')
    plt.plot(aas,delta4s,'m',label='Newtonian RD')
    plt.plot(aas,delta2s,'r',label='MD')
    plt.plot(aas,delta3s,'g',label='Synchronous RD')
    plt.plot([gt.aeq,gt.aeq],[np.min([delta1s,delta2s,delta3s]),np.max([delta1s,delta2s,delta3s])],'c--',label=r'$a_{\mathrm{eq}}$')
    plt.plot([gt.a(zinit),gt.a(zinit)],[np.min([delta1s,delta2s,delta3s]),np.max([delta1s,delta2s,delta3s])],'m--',label=r'initial condition')
    plt.plot([gt.a(zH),gt.a(zH)],[np.min([delta1s,delta2s,delta3s]),np.max([delta1s,delta2s,delta3s])],'y--',label=r'horizon entry')
    plt.title(r'$k = %s$ %s'%(latex_float(kplot),gt.unitTexMpc1))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$a$')
    plt.ylabel(r'$\delta$')
    plt.ylim(delta4s.min()/2,delta1s.max())
    plt.xlim(aas.min(),aas.max())
    #plt.ylim(-.1,.3)
    plt.legend(loc=2)
      
  if plot_hor:
    zH = (-gt.OmegaM+np.sqrt(gt.OmegaM**2+4.*gt.OmegaR*(kplot/gt.H0_in_Mpc)**2))/(2.*gt.OmegaR)-1.
    aas = np.logspace(np.log10(gt.a(zH))-.2,-1,200)
    horizons = gt.horizon(gt.z(aas))
    phorizons = np.zeros(aas.size)
    for i in range(aas.size):
      phorizons[i] = gt.particleHorizon(gt.z(aas[i]))
    kzeros = gt.kzero(aas)
    plt.figure(6)
    plt.clf()
    plt.plot(aas,horizons,label=r'horizon $1/aH$')
    plt.plot(aas,phorizons,label=r'particle horizon')
    plt.plot(aas,1./kzeros,label='Meszaros soln breakdown')
    plt.plot(aas,np.full(aas.size,1./kplot),'k-',label=r'$1/k_\mathrm{halo}$')
    plt.plot([gt.aeq,gt.aeq],[1e-4*np.min(horizons),1e4*np.max(horizons)],'c--',label=r'$a_{\mathrm{eq}}$')
    plt.plot([gt.a(zinit),gt.a(zinit)],[1e-4*np.min(horizons),1e4*np.max(horizons)],'m--',label=r'initial condition')
    plt.plot([gt.a(zH),gt.a(zH)],[1e-4*np.min(horizons),1e4*np.max(horizons)],'y--',label=r'horizon entry')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$a$')
    plt.ylabel(r'length (%s)'%gt.unitTexMpc)
    plt.ylim(horizons.min()/10,phorizons.max())
    plt.legend(loc=2)
      
  if plot_prim:
    plt.figure(10,figsize=(5,2.5))
    plt.clf()
    PR = gt.R0(k)**2
    #plt.loglog(k,PR,'k--')
    #plt.loglog(k,PR,'k--',label='')
    modify(k,PR)
    plt.loglog(k,PR,'k',label='primordial curvature')
    #plt.legend(loc=2)
    plt.loglog(np.full(2,kminBox),[1e-100,1e100],'k--')
    plt.loglog(np.full(2,kmaxReq),[1e-100,1e100],'k--')
    plt.ylim(gt.A_s/6,PR.max()*3)
    plt.xlabel(r'$k$ (%s)'%gt.unitTexMpc1)
    plt.ylabel(r'$\mathcal{P}_\zeta(k)$')
    plt.tight_layout()
  
  plt.show()
#!/usr/bin/env python
import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
import os
#import copy as cp
#import time
import sys
#from scipy.optimize import minimize
#import random
#import scipy.stats as st
#import xlsxwriter

#path where you put func_ORCHIPRIMv0.py file
#if the func_ORCHIPRIMv0.y is in the same directory of this script, you use line in the following line
func_path=os.getcwd()
pydir=func_path
sys.path.append(pydir)
from ORCHIMIC import microbe_growth

#starttime=time.time()

par_path=func_path+'/par/'

#path to save output
outdir0='/home/surface3/yhuang/microbe/'

#some flags to check the model, for example whether it generate unreasonable values like native value or nan or inf, and mass blance in the model
OK_test='y' #A test run, save the output into your test path
OK_check='n'  #check whether the model will generate some unreasonable values
if OK_check=='y':
  OK_print='y'
else:
  OK_print='n'

#ts=24 #unit hour
OK_constant_CAE='y' #not used in this verion, but give chance to use CAE that decrease with temeprature

ts=24  #time step of the model, default 24
OK_N='n'  #y:with N dynamics; n: without N dynamics
RLRS_method=1 #default:1 (used in paper); 1: based on AvailC pool; 2: use new decomposed C first
cost_method=1 #defaul:1 (used in paper)

vers='v0' #version number
model=[1]  #define the ORCHIMIC variants you want to use: 1: one general MFT; 2: one FOM specialist + one SOM specialist; 3: one FOM specialist + one SOM specialist + one general MFT

# default is 1, should not be changed
global N
N=1
#molar mass for C, N and P
MMC=12.011
MMN=14.007
MMP=30.973761998

#default setting, set C/N ratio for different microbial function types (MFTs)
BCN=np.zeros((4,1))
BCN[:,0]=np.array([6.12,4.59,8.3,6.12]) #Mouginot et al., 2014
ECN=BCN

#recyle fraction of dead microbes to Avail pool
sC=0.06 #0.06 from C. Kaiser et al. 2014 and 2015
sN=0.06
sP=0.06

OK_vegNuptake='n' #activating uptake of N by plants or not
OK_loss='n'       #activating leaching or not
OK_adsorb_C='y'   #activating desorption and adsorption for C or not
OK_adsorb_NP='y'  #activating desorption and adsorption for N (P will added in future version) or not
OK_control='n'    #if yes, the modeled the concentration of microbes will be reset to a minimum value if the modeled the concentration of microbes is lower than the prescribed minimum value
OK_restart='n'    #restart run or not 
OK_reset_Lr='n'   #reset FOM derived C fraction in each pool or not
exp='test'        #name of experiment
xxnv=22           #number of parameter values read from other file

#fixed parameter values got from CENTURY model
SStoSP=0.03       #Flux go to SP from SS due to physichemical protection
SAtoSP=0.004      #Flux go to SP from SA due to physichemical protection

#Character of FOM
LCin=1.6e-4 #input rate, 1.6e-4 mg C/g soil/h from Wang et al., 2013
lignin=0.2+np.zeros((N)) #lignin content of FOM
LCNin=50.+np.zeros((N)) #C/N ratio of FOM, McGroddy et al., 2004 

#Temperature, constant temperature during running
T=273.15+22

for MFT in model:
    if OK_test=='y':
      N=1
      outdir=outdir0+'tranp/test/'
    else:
      outdir=outdir0+'tranp/'+vers+'/T'+str(T)+'/MFT'+str(MFT)+'/'+exp+'/'
    
    if not os.path.exists(outdir):
      os.makedirs(outdir)
    
    #if it is a restart simulation,provide the path of output for last year
    #if OK_restart=='y':
      #restart path
      #restart year
      #restart_dir=restart_dir
    print 'use model MFT '+str(MFT)
    print 'N=',N
    
    if MFT==1:
            seedn=0
            mn=1
            MFTtype=[3] #[1]*3+[2]*3+[3]*3+[4]*3+[5]*3+[6]*3
    elif MFT==2:
            seedn=1
            mn=2
            MFTtype=[1,2]
    elif MFT==3:
            seedn=3
            mn=3
            MFTtype=[3,1,2]

    if OK_N=='y':
      xt=np.load(par_path+'/CN-MFT'+str(MFT)+'.npy')
    elif OK_N=='n':
      xt=np.load(par_path+'/C-MFT'+str(MFT)+'.npy')

    xx=np.zeros((np.size(xt),1))
    xx[:,0]=xt
    
    #some values that need to default
    H=0.5
    pH=7.
    Actf=0.05
    Pasf=0.5
    clay=0.2
    
    
    #Default values for parameters before optimization, some of them may not be used, because they are read from files
    Adj_GL          =1000.
    Adj_LS          =1.        #5.02/0.96
    Adj_SA          =37.
    Adj_SP          =29.	
    AvailCrin       =0.098/24. #0.002 #0.005773157
    b               =0.001
    BArin           =0.25+np.zeros((N))
    CAE             =0.6
    dENZ            =0.001
    dMFTref         =0.001
    Ea_main         =20.       #60. #Tang and Riley, 2015 #20.     #see Wang et al., 2013
    EaGL            =28.9      #Odebunmi and Owalude, 2007
    EaKM            =30.       #Davidson et al. (2006) from Wang et al., 2013
    EaLM            =37.       #Wang et al., 2012b from Wang et al., 2013
    EaLS            =53.       #Wang et al., 2012b from Wang et al., 2013
    Ea_mob	    =20.       #Check references cited by Wang et al., 2013
    EaSA            =47.
    EaSP            =47.
    EaSS            =47.       #Allison et al., 2010
    Ea_stab	    =5.        #Check references cited by Wang et al., 2013
    Ea_uptake       =47.       #47 from Allison et al., 2010
    expc            =1.+np.zeros((N))
    efLS            =1.+np.zeros((N))
    efSA            =1.+np.zeros((N))
    efSP            =1.+np.zeros((N))
    efSS            =1.+np.zeros((N))
    ELrin           =np.array([0.00,0.25,0.75,0.50]) #fraction of enzyme production capacity can be potentially used to producing FOM decomposing enzymes
    FErin	    =0.05
    H0              =H         #0.6
    Hs              =0.4
    KBAC            =6.        #6(1-11) g soil/mg C Wang et al., 2013
    K_CAE           =0.012     #Devevre and Horwath 2000 from Wang et al., 2013
    KdesC           =0.001+np.zeros(N)   #Wang et al., 2013
    Kdin            =0.
    Kein            =1e-4
    Ke_min          =0.1+np.zeros((N))
    KlossC          =0.05/365./24.+np.zeros((N))
    KlossN          =0.05/365/24+np.zeros((N))
    KlossN2         =0.05+np.zeros((N)) #KlossN2 is totally different from KlossC, as in Parton et al., 1987, %5 of totoal mineralization flux as volatilization
    KME_Avail       =0.26
    KME_FOM         =50.
    KME_SOM         =250.
    Kr_ref          =0.01
    KMLC            =50.
    KMSC            =250.
    KM_uptake       =0.26
    LMtoSS          =0.+np.zeros((N))
    LStoSS          =0.+np.zeros((N)) 
    pH0             =pH          #pH related parameters cannot be optimized according current data we used in the paper
    pH0_dec         =pH      
    pHs		    =2.
    pHs_dec         =2.
    SACrin	    =Actf
    sC              =0.06
    SErin	    =0.05
    sN              =0.06
    sP              =0.06
    SPCrin	    =Pasf
    AdsorbCmax	    =1.7         #mmol C/kg soil from 1.7 mg C/g soil from Mayes et al. (2012) in wang et al., 2013   
    T0              =T
    Ts		    =999.
    VmaxLMC         =2.5
    VmaxSSC         =1.
    VmaxUptake_refin=0.24
    
    
    ##Read values from outside of function for parameters           
    ll=-1
    #ll=ll+1; Adj_GL          =xx[ll]
    #ll=ll+1; Adj_LS          =xx[ll]
    ll=ll+1; Adj_SA          =xx[ll]
    ll=ll+1; Adj_SP          =xx[ll]
    ll=ll+1; AvailCrin       =xx[ll]
    ll=ll+1; b               =xx[ll]
    ll=ll+1; BArin           =xx[ll]
    ll=ll+1; CAE             =xx[ll]
    #ll=ll+1; dENZ            =xx[ll]
    ll=ll+1; dMFTref         =xx[ll]
    #ll=ll+1; Ea_main         =xx[ll]
    #ll=ll+1; EaGL            =xx[ll]
    #ll=ll+1; EaLM            =xx[ll]
    #ll=ll+1; EaLS            =xx[ll]
    #ll=ll+1; EaSA            =xx[ll]
    #ll=ll+1; EaSP            =xx[ll]
    #ll=ll+1; EaSS            =xx[ll]
    #ll=ll+1; Ea_uptake       =xx[ll]
    ll=ll+1; expc           =xx[ll]
    #ll=ll+1; efLS           =xx[ll]
    #ll=ll+1; efSA           =xx[ll]
    #ll=ll+1; efSP           =xx[ll]
    #ll=ll+1; efSS           =xx[ll]
    ll=ll+1; FErin		=xx[ll]
    #ll=ll+1; H0              =xx[ll]
    #ll=ll+1; Hs              =xx[ll]
    #ll=ll+1;KabsC		=xx[ll]
    ll=ll+1; KBAC            =xx[ll]
    #ll=ll+1; K_CAE           =xx[ll]
    #ll=ll+1; Kdin           =xx[ll]
    ll=ll+1; KdesC           =xx[ll]
    ll=ll+1; Kein            =xx[ll]
    #ll=ll+1; Ke_min          =xx[ll]
    #ll=ll+1; KlossC          =xx[ll]
    #ll=ll+1; KME_Avail      =xx[ll]
    #ll=ll+1; KME_FOM        =xx[ll]
    #ll=ll+1; KME_SOM        =xx[ll]
    ll=ll+1; Kr_ref          =xx[ll]
    ll=ll+1; KMLC            =xx[ll]
    ll=ll+1; KMSC            =xx[ll]
    ll=ll+1; KM_uptake       =xx[ll]
    ll=ll+1; LStoSS         =xx[ll]
    #ll=ll+1; pH0_dec         =xx[ll]
    #ll=ll+1; pHs_dec         =xx[ll]
    #ll=ll+1; SACrin		=xx[ll]
    ll=ll+1; SErin          =xx[ll]
    #ll=ll+1; SPCrin		 =xx[ll]
    ll=ll+1; AdsorbCmax	=xx[ll]
    #ll=ll+1; T0              =xx[ll]
    #ll=ll+1; Ts              =xx[ll]
    ll=ll+1; VmaxLMC         =xx[ll]
    ll=ll+1; VmaxSSC         =xx[ll]
    ll=ll+1; VmaxUptake_refin=xx[ll]
    
    VmaxLMC=VmaxLMC/Adj_LS #/2.
    Adj_LS=1.
    
    KabsC=KdesC*KBAC
    
    KME_FOM=KMLC
    KME_SOM=KMSC
    KME_Avail=KM_uptake
    
    if not (ll+1)==xxnv:
      print 'Error: read values for parameters'
      print 'll+1=',ll+1, '  xxnv=',xxnv
      sys.exit()
    
    #set initial valuse 
    TC=20.    #initial total soil cabron stock
    LC0=1.    #initial litter carbon stock
    SCN=12.   #initial soil organic matter C/N ratio
    
    #ANin=ANin*1000./MMN
    #APin=APin*1000./MMP
    
    AvailC0=TC*0.005+np.zeros((1,N))
    #AdsorbC0=TC*AdsorbCrin+np.zeros(N)
    AdsorbC0=AvailC0*0.
    
    AvailN0=AvailC0/SCN
    AdsorbN0=AdsorbC0/SCN
    
    #Initialization of pools
    BA0=TC*0.025/mn*0.3+np.zeros((mn,N))
    BD0=TC*0.025/mn*0.7+np.zeros((mn,N))
    
    ELC0=BA0*0.1
    ESC0=BA0*0.1
    ##ELC0=(BA0+BD0)*FEr_max*FErin
    ##ESC0=(BA0+BD0)*SEr_max*SErin
    
    for m in np.arange(mn):
      if MFTtype[m]==0: #cheater as a MFT, just for test 
        ELC0[m,:]=0.
        ESC0[m,:]=0.
    
    SAC0=TC*0.05+np.zeros((1,N))
    SPC0=TC*0.5+np.zeros((1,N))
    SSC0=TC*0.45+np.zeros((1,N))
    
    SAN0=SAC0/SCN
    SPN0=SPC0/SCN
    SSN0=SSC0/SCN
    
    #AvailCin=np.zeros((1,N))
    #AvailCin=ACin+np.zeros(N)
    #AvailNin=ANin+np.zeros(N)
    #AvailPin=APin+np.zeros(N)
    GL0=0.+np.zeros((1,N)) #ACin+np.zeros((1,N))
    #LC0=LCin+np.zeros((1,N))
    
    AvailCLr0=0.+np.zeros((1,N))
    #AvailCLr0=ACin/(ACin+AvailC0)+np.zeros((1,N))
    AvailCSr0=1.-AvailCLr0
    AdsorbCLr0=0.+np.zeros((1,N))
    AdsorbCSr0=1.-AdsorbCLr0
    SACLr0=0.+np.zeros((1,N))
    SACSr0=1.-SACLr0
    SSCLr0=0.+np.zeros((1,N))
    SSCSr0=1.-SSCLr0
    SPCLr0=0.+np.zeros((1,N))
    SPCSr0=1.-SPCLr0
    BALr0=0.+np.zeros((mn,N))
    BASr0=1.- BALr0
    BDLr0=0.+np.zeros((mn,N))
    BDSr0=1.- BDLr0
    ELCLr0=0.+np.zeros((mn,N))
    ELCSr0=1.- ELCLr0
    #ELNLr0=0.+np.zeros((mn,N))
    #ELNSr0=1.- ELNLr0
    #ELPLr0=0.+np.zeros((mn,N))
    #ELPSr0=1.- ELPLr0
    ESCLr0=0.+np.zeros((mn,N))
    ESCSr0=1.- ESCLr0
    #ESNLr0=0.+np.zeros((mn,N))
    #ESNSr0=1.- ESNLr0
    #ESPLr0=0.+np.zeros((mn,N))
    #ESPSr0=1.- ESPLr0
    
    LMC0=LC0*0.5+np.zeros((1,N)) #*LM_fraction
    #LitterMN0=LitterN0*(1.-lignin);
    #LitterMP0=LitterP0*(1.-lignin);
    LSC0=LC0*0.5+np.zeros((1,N)) #*(1-LM_fraction);
    LLf0=LSC0*0.+0.5 #LC0*lignin/LSC0
    LLf0[LSC0==0]=0.
    
    LMN0=LMC0/50.
    LSN0=LSC0/150.
    
    if OK_restart=='y':
      GLtmp=np.load(restart_dir+'GL.npy')
      LMCtmp=np.load(restart_dir+'LMC.npy')
      LSCtmp=np.load(restart_dir+'LSC.npy')
      LMNtmp=np.load(restart_dir+'LMN.npy')
      LSNtmp=np.load(restart_dir+'LSN.npy')
      LLftmp=np.load(restart_dir+'LLf.npy')
      
      SACtmp=np.load(restart_dir+'SAC.npy')
      SSCtmp=np.load(restart_dir+'SSC.npy')
      SPCtmp=np.load(restart_dir+'SPC.npy')
      SANtmp=np.load(restart_dir+'SAN.npy')
      SSNtmp=np.load(restart_dir+'SSN.npy')
      SPNtmp=np.load(restart_dir+'SPN.npy')
    
      AvailCtmp=np.load(restart_dir+'AvailC.npy')
      AdsorbCtmp=np.load(restart_dir+'AdsorbC.npy')
      AvailNtmp=np.load(restart_dir+'AvailN.npy')
      AdsorbNtmp=np.load(restart_dir+'AdsorbN.npy')
    
      BAtmp=np.load(restart_dir+'BA.npy')
      BDtmp=np.load(restart_dir+'BD.npy')
      ELCtmp=np.load(restart_dir+'ELC.npy')
      ESCtmp=np.load(restart_dir+'ESC.npy')
    
      GL0[0,:]=GLtmp[-1,:]
      LMC0[0,:]=LMCtmp[-1,:]
      LSC0[0,:]=LSCtmp[-1,:]
      LMN0[0,:]=LMNtmp[-1,:]
      LSN0[0,:]=LSNtmp[-1,:]
      LLf0[0,:]=LLftmp[-1,:]
    
      SAC0[0,:]=SACtmp[-1,:]
      SSC0[0,:]=SSCtmp[-1,:]
      SPC0[0,:]=SPCtmp[-1,:]
      SAN0[0,:]=SANtmp[-1,:]
      SSN0[0,:]=SSNtmp[-1,:]
      SPN0[0,:]=SPNtmp[-1,:]
    
      AvailC0[0,:]=AvailCtmp[-1,:]
      AdsorbC0[0,:]=AdsorbCtmp[-1,:]
      AvailN0[0,:]=AvailNtmp[-1,:]
      AdsorbN0[0,:]=AdsorbNtmp[-1,:]
    
      BA0[:,:]=BAtmp[:,-1,:]
      BD0[:,:]=BDtmp[:,-1,:]
      ELC0[:,:]=ELCtmp[:,-1,:]
      ESC0[:,:]=ESCtmp[:,-1,:]
    
    
    yy1=0
    yy2=1000
    dd0=365
    #LCin=1.6e-4 #1.6e-4 mg C/g soil/h from Wang et al., 2013
    #lignin=0.2+np.zeros((N))
    #LCNin=50.+np.zeros((N)) #McGroddy et al., 2004  
    ACin=np.zeros((N))
    ANin=np.zeros((N))
    for yy in range(yy1,yy2,1):
        #print 'yy=',yy
        
        #if yy==10 and not exp=='spinup':
        if yy==10:
          if exp=='5K' or exp=='5K_2I':
            T=T+5.
          if exp=='2I' or exp=='5K_2I':
            LCin=LCin*2.
        T0=T
        outdirtmp=outdir+str(yy)+'/'
        if not os.path.exists(outdirtmp):
            os.makedirs(outdirtmp)
        t0=dd0*24/ts
        #initialization
        GL=np.zeros((t0+1,N))
        #VEGC=np.zeros((t0+1,N))
        #VEGN=np.zeros((t0+1,N))
        #VEGP=np.zeros((t0+1,N))
        #VEGCg=np.zeros((t0,N))
        #VEGNg=np.zeros((t0,N))
        #VEGPg=np.zeros((t0,N))
        BA=np.zeros((mn,t0+1,N));
        BD=np.zeros((mn,t0+1,N));
        MCUE1=np.zeros((t0,N));
        MCUE2=np.zeros((t0,N));
        ELC=np.zeros((mn,t0+1,N));
        ESC=np.zeros((mn,t0+1,N));
        #ELN=np.zeros((mn,t0+1,N));
        #ESN=np.zeros((mn,t0+1,N));
        #ELP=np.zeros((mn,t0+1,N));
        #ESP=np.zeros((mn,t0+1,N));
        ELCeq=np.zeros((t0,N));
        ESCeq=np.zeros((t0,N));
        #ELNeq=np.zeros((t0,N));
        #ESNeq=np.zeros((t0,N));
        #ELPeq=np.zeros((t0,N));
        #ESPeq=np.zeros((t0,N));
        LMC=np.zeros((t0+1,N));
        LSC=np.zeros((t0+1,N));
        LLf=np.zeros((t0+1,N));
        SAC=np.zeros((t0+1,N));
        SSC=np.zeros((t0+1,N));
        SPC=np.zeros((t0+1,N));
        AvailC=np.zeros((t0+1,N));
        AdsorbC=np.zeros((t0+1,N));
        
        LMN=np.zeros((t0+1,N));
        LSN=np.zeros((t0+1,N));
        SAN=np.zeros((t0+1,N));
        SSN=np.zeros((t0+1,N));
        SPN=np.zeros((t0+1,N));
        AvailN=np.zeros((t0+1,N));
        AdsorbN=np.zeros((t0+1,N));
    
        saturation_ratio=np.zeros((mn,t0,N)); #
        gC=np.zeros((mn,t0,N)); #growth rate at each time step
        gN=np.zeros((mn,t0,N)); #growth rate at each time step
        #gP=np.zeros((mn,t0,N));
        g=np.zeros((mn,t0,N));
        Rm=np.zeros((mn,t0,N)); #maintenance respiration at each time step
        Rg=np.zeros((mn,t0,N)); #growth respiration at each time step
        Rd=np.zeros((t0,N));
        Ro=np.zeros((mn,t0,N));
        RC=np.zeros((t0,N));
        RL=np.zeros((t0,N));
        RS=np.zeros((t0,N));
    
        dGL=np.zeros((t0,N));
        dLMC=np.zeros((t0,N));
        #dLMN=np.zeros((t0,N));
        #dLMP=np.zeros((t0,N));
        dLSC=np.zeros((t0,N));
        #dLSN=np.zeros((t0,N));
        #dLSP=np.zeros((t0,N));
        dSAC=np.zeros((t0,N));
        #dSAN=np.zeros((t0,N));
        #dSAP=np.zeros((t0,N));
        dSSC=np.zeros((t0,N));
        #dSSN=np.zeros((t0,N));
        #dSSP=np.zeros((t0,N));
        dSPC=np.zeros((t0,N));
        #dSPN=np.zeros((t0,N));
        #dSPP=np.zeros((t0,N));
        BAg=np.zeros((mn,t0,N));
        BAd=np.zeros((mn,t0,N));
        Rma=np.zeros((mn,t0,N));
        Rmd=np.zeros((mn,t0,N));
        AvailCstab=np.zeros((t0,N));
        AvailCmob=np.zeros((t0,N));
        AvailCloss=np.zeros((t0,N));
        #AvailNstab=np.zeros((t0,N));
        #AvailNmob=np.zeros((t0,N));
        #AvailNloss=np.zeros((t0,N));
        #AvailPstab=np.zeros((t0,N));
        #AvailPmob=np.zeros((t0,N));
        #AvailPloss=np.zeros((t0,N));
        ELCg=np.zeros((mn,t0,N));
        #ELNg=np.zeros((mn,t0,N));
        #ELPg=np.zeros((mn,t0,N));
        ESCg=np.zeros((mn,t0,N));
        #ESNg=np.zeros((mn,t0,N));
        #ESPg=np.zeros((mn,t0,N));
    
        AvailCLr=np.zeros((t0+1,N));
        AvailCSr=np.zeros((t0+1,N));
        AdsorbCLr=np.zeros((t0+1,N));
        AdsorbCSr=np.zeros((t0+1,N));
        SACLr=np.zeros((t0+1,N));
        SACSr=np.zeros((t0+1,N));
        SSCLr=np.zeros((t0+1,N));
        SSCSr=np.zeros((t0+1,N));
        SPCLr=np.zeros((t0+1,N));
        SPCSr=np.zeros((t0+1,N));
        BALr=np.zeros((mn,t0+1,N))
        BASr=np.zeros((mn,t0+1,N))
        BDLr=np.zeros((mn,t0+1,N))
        BDSr=np.zeros((mn,t0+1,N))
        ELCLr=np.zeros((mn,t0+1,N))
        ELCSr=np.zeros((mn,t0+1,N))
        #ELNLr=np.zeros((mn,t0+1,N))
        #ELNSr=np.zeros((mn,t0+1,N))
        #ELPLr=np.zeros((mn,t0+1,N))
        #ELPSr=np.zeros((mn,t0+1,N))
        ESCLr=np.zeros((mn,t0+1,N))
        ESCSr=np.zeros((mn,t0+1,N))
        #ESNLr=np.zeros((mn,t0+1,N))
        #ESNSr=np.zeros((mn,t0+1,N))
        #ESPLr=np.zeros((mn,t0+1,N))
        #ESPSr=np.zeros((mn,t0+1,N))
        
        if yy==yy1:
            BA[:,0]=BA0;
            BD[:,0]=BD0;
            ELC[:,0]=ELC0;
            ESC[:,0]=ESC0;
            #ELN[:,0]=ELN0;
            #ESN[:,0]=ESN0;
            #ELP[:,0]=ELP0;
            #ESP[:,0]=ESP0;
            GL[0]=GL0
            LMC[0]=LMC0;
            LSC[0]=LSC0;
            LLf[0]=LLf0
            SAC[0]=SAC0;
            SSC[0]=SSC0;
            SPC[0]=SPC0;
            AvailC[0]=AvailC0;
            AdsorbC[0]=AdsorbC0
         
            LMN[0]=LMN0;
            LSN[0]=LSN0;
            SAN[0]=SAN0;
            SSN[0]=SSN0;
            SPN[0]=SPN0;
            AvailN[0]=AvailN0;
            AdsorbN[0]=AdsorbN0
    
            AvailCLr[0]=AvailCLr0;
            AvailCSr[0]=AvailCSr0;
            AdsorbCLr[0]=AdsorbCLr0;
            AdsorbCSr[0]=AdsorbCSr0;
            SACLr[0]=SACLr0
            SACSr[0]=SACSr0
            SSCLr[0]=SSCLr0
            SSCSr[0]=SSCSr0
            SPCLr[0]=SPCLr0
            SPCSr[0]=SPCSr0
            BALr[:,0]=BALr0
            BASr[:,0]=BASr0
            BDLr[:,0]=BDLr0
            BDSr[:,0]=BDSr0
            ELCLr[:,0]=ELCLr0
            ELCSr[:,0]=ELCSr0
            #ELNLr[:,0]=ELNLr0
            #ELNSr[:,0]=ELNSr0
            #ELPLr[:,0]=ELPLr0
            #ELPSr[:,0]=ELPSr0
            ESCLr[:,0]=ESCLr0
            ESCSr[:,0]=ESCSr0
            #ESNLr[:,0]=ESNLr0
            #ESNSr[:,0]=ESNSr0
            #ESPLr[:,0]=ESPLr0
            #ESPSr[:,0]=ESPSr0
        else:
            BA[:,0]=BA_out;
            BD[:,0]=BD_out;
            ELC[:,0]=ELC_out;
            ESC[:,0]=ESC_out;
            GL[0]=GL_out
            LMC[0]=LMC_out;
            LSC[0]=LSC_out
            LLf[0]=LLf_out
            SAC[0]=SAC_out
            SSC[0]=SSC_out
            SPC[0]=SPC_out
            AvailC[0]=AvailC_out
            AdsorbC[0]=AdsorbC_out
            
            LMN[0]=LMN_out;
            LSN[0]=LSN_out
            SAN[0]=SAN_out
            SSN[0]=SSN_out
            SPN[0]=SPN_out
            AvailN[0]=AvailN_out
            AdsorbN[0]=AdsorbN_out
    
            ff=1.
            if OK_reset_Lr=='y':
              ff=0.
            
            AvailCLr[0]=AvailCLr_out*ff
            AvailCSr[0]=AvailCSr_out*ff
            AdsorbCLr[0]=AdsorbCLr_out*ff
            AdsorbCSr[0]=AdsorbCSr_out*ff
            SACLr[0]=SACLr_out*ff
            SACSr[0]=SACSr_out*ff
            SSCLr[0]=SSCLr_out*ff
            SSCSr[0]=SSCSr_out*ff
            SPCLr[0]=SPCLr_out*ff
            SPCSr[0]=SPCSr_out*ff
            BALr[:,0]=BALr_out*ff
            BASr[:,0]=BASr_out*ff
            BDLr[:,0]=BDLr_out*ff
            BDSr[:,0]=BDSr_out*ff
            ELCLr[:,0]=ELCLr_out*ff
            ELCSr[:,0]=ELCSr_out*ff
            ESCLr[:,0]=ESCLr_out*ff
            ESCSr[:,0]=ESCSr_out*ff
    
        tt=0
        for dd in range(dd0):
          for t in range(24/ts):
            print 'yy=',yy,'dd=',dd,'t=',t
            #LCin=1.6e-4 #1.6e-4 mg C/g soil/h from Wang et al., 2013
            #lignin=0.2+np.zeros((N))
            #LCNin=50.+np.zeros((N)) #McGroddy et al., 2004  
            #ACin=np.zeros((N))  
            #ANin=np.zeros((N))
            if yy==yy1 and tt==0:
                GL_in=GL0
                LMC_in=LMC0
                LSC_in=LSC0
                LLf_in=LLf0
                LMN_in=LMN0
                LSN_in=LSN0
                SAC_in=SAC0
                SSC_in=SSC0
                SPC_in=SPC0
                SAN_in=SAN0
                SSN_in=SSN0
                SPN_in=SPN0
                AvailC_in=AvailC0
                AdsorbC_in=AdsorbC0
                AvailN_in=AvailN0
                AdsorbN_in=AdsorbN0
                BA_in=BA0
                BD_in=BD0
                ELC_in=ELC0
                ESC_in=ESC0
                AvailCLr_in=AvailCLr0
                AvailCSr_in=AvailCSr0
                AdsorbCLr_in=AdsorbCLr0
                AdsorbCSr_in=AdsorbCSr0
                SACLr_in=SACLr0
                SACSr_in=SACSr0
                SSCLr_in=SSCLr0
                SSCSr_in=SSCSr0
                SPCLr_in=SPCLr0
                SPCSr_in=SPCSr0
                BALr_in=BALr0
                BASr_in=BASr0
                BDLr_in=BDLr0
                BDSr_in=BDSr0
                ELCLr_in=ELCLr0
                ELCSr_in=ELCSr0
                ESCLr_in=ESCLr0
                ESCSr_in=ESCSr0
            else:
                GL_in=np.array(GL_out)
                LMC_in=np.array(LMC_out)
                LSC_in=np.array(LSC_out)
                LLf_in=np.array(LLf_out)
                LMN_in=np.array(LMN_out)
                LSN_in=np.array(LSN_out)
                SAC_in=np.array(SAC_out)
                SSC_in=np.array(SSC_out)
                SPC_in=np.array(SPC_out)
                SAN_in=np.array(SAN_out)
                SSN_in=np.array(SSN_out)
                SPN_in=np.array(SPN_out)
                AvailC_in=np.array(AvailC_out)
                AdsorbC_in=np.array(AdsorbC_out)
                AvailN_in=np.array(AvailN_out)
                AdsorbN_in=np.array(AdsorbN_out)
                BA_in=np.array(BA_out)
                BD_in=np.array(BD_out)
                ELC_in=np.array(ELC_out)
                ESC_in=np.array(ESC_out)
                AvailCLr_in=np.array(AvailCLr_out)
                AvailCSr_in=np.array(AvailCSr_out)
                AdsorbCLr_in=np.array(AdsorbCLr_out)
                AdsorbCSr_in=np.array(AdsorbCSr_out)
                SACLr_in=np.array(SACLr_out)
                SACSr_in=np.array(SACSr_out)
                SSCLr_in=np.array(SSCLr_out)
                SSCSr_in=np.array(SSCSr_out)
                SPCLr_in=np.array(SPCLr_out)
                SPCSr_in=np.array(SPCSr_out)
                BALr_in=np.array(BALr_out)
                BASr_in=np.array(BASr_out)
                BDLr_in=np.array(BDLr_out)
                BDSr_in=np.array(BDSr_out)
                ELCLr_in=np.array(ELCLr_out)
                ELCSr_in=np.array(ELCSr_out)
                ESCLr_in=np.array(ESCLr_out)
                ESCSr_in=np.array(ESCSr_out)
            
            GL_out,LMC_out,LMN_out,LSC_out,LSN_out,LLf_out,SAC_out,SAN_out,SSC_out,SSN_out,SPC_out,SPN_out,AvailC_out,AvailN_out,BA_out,BD_out,ELC_out,ESC_out,ELCeq_out,ESCeq_out,AdsorbC_out,AdsorbN_out,MCUE1_out,MCUE2_out,saturation_ratio_out,Rd_out,Ro_out,Rg_out,Rm_out,RC_out,gC_out,gN_out,g_out,AvailCLr_out,AvailCSr_out,AdsorbCLr_out,AdsorbCSr_out,SACLr_out,SACSr_out,SSCLr_out,SSCSr_out,SPCLr_out,SPCSr_out,BALr_out,BASr_out,BDLr_out,BDSr_out,ELCLr_out,ELCSr_out,ESCLr_out,ESCSr_out,dGL_out,dLMC_out,dLMN_out,dLSC_out,dLSN_out,dSAC_out,dSAN_out,dSSC_out,dSSN_out,dSPC_out,dSPN_out,BAg_out,BAd_out,Rma_out,Rmd_out,RL_out,RS_out,AvailCstab_out,AvailNstab_out,AvailCmob_out,AvialNmob_out,AvailCloss_out,AvailNloss_out,ELCg_out,ESCg_out=microbe_growth(OK_check,OK_print,OK_loss,OK_N,OK_vegNuptake,N,ts,RLRS_method,K_CAE,VmaxUptake_refin,Kdin,dMFTref,Kein,Kr_ref,b,dENZ,VmaxLMC,KMLC,Ea_stab,Ea_mob,EaGL,EaKM,EaLM,EaLS,VmaxSSC,KMSC,EaSS,AdsorbCmax,KdesC,KabsC,CAE,Adj_GL,Adj_LS,Adj_SA,Adj_SP,Ke_min,Ea_uptake,Ea_main,EaSA,EaSP,KM_uptake,ELrin,pH0,T0,H0,pHs,Ts,Hs,pH0_dec,pHs_dec,KlossC,KlossN,KlossN2,LCin,lignin,LCNin,ACin,ANin,GL_in,LMC_in,LMN_in,LSC_in,LSN_in,LLf_in,SAC_in,SAN_in,SSC_in,SSN_in,SPC_in,SPN_in,AvailC_in,AvailN_in,AdsorbC_in,AdsorbN_in,BA_in,BD_in,MFTtype,BCN,ELC_in,ESC_in,ECN,OK_control,OK_constant_CAE,OK_adsorb_C,T,H,pH,clay,mn,sC,sN,AvailCLr_in,AvailCSr_in,AdsorbCLr_in,AdsorbCSr_in,SACLr_in,SACSr_in,SSCLr_in,SSCSr_in,SPCLr_in,SPCSr_in,BALr_in,BASr_in,BDLr_in,BDSr_in,ELCLr_in,ELCSr_in,ESCLr_in,ESCSr_in,LMtoSS,LStoSS,SStoSP,SAtoSP,expc,efLS,efSA,efSS,efSP,KME_FOM,KME_SOM,KME_Avail)
            
            BA[:,tt+1]=BA_out;
            BD[:,tt+1]=BD_out;
            ELC[:,tt+1]=ELC_out;
            ESC[:,tt+1]=ESC_out;
            GL[tt+1]=GL_out
            LMC[tt+1]=LMC_out;
            LSC[tt+1]=LSC_out
            LLf[tt+1]=LLf_out
            SAC[tt+1]=SAC_out
            SSC[tt+1]=SSC_out
            SPC[tt+1]=SPC_out
            AvailC[tt+1]=AvailC_out
            AdsorbC[tt+1]=AdsorbC_out
            
            LMN[tt+1]=LMN_out;
            LSN[tt+1]=LSN_out
            SAN[tt+1]=SAN_out
            SSN[tt+1]=SSN_out
            SPN[tt+1]=SPN_out
            AvailN[tt+1]=AvailN_out
            AdsorbN[tt+1]=AdsorbN_out
            
            AvailCLr[tt+1]=AvailCLr_out
            AvailCSr[tt+1]=AvailCSr_out
            AdsorbCLr[tt+1]=AdsorbCLr_out
            AdsorbCSr[tt+1]=AdsorbCSr_out;
            SACLr[tt+1]=SACLr_out
            SACSr[tt+1]=SACSr_out
            SSCLr[tt+1]=SSCLr_out
            SSCSr[tt+1]=SSCSr_out
            SPCLr[tt+1]=SPCLr_out
            SPCSr[tt+1]=SPCSr_out
            BALr[:,tt+1]=BALr_out
            BASr[:,tt+1]=BASr_out
            BDLr[:,tt+1]=BDLr_out
            BDSr[:,tt+1]=BDSr_out
            ELCLr[:,tt+1]=ELCLr_out
            ELCSr[:,tt+1]=ELCSr_out
            ESCLr[:,tt+1]=ESCLr_out
            ESCSr[:,tt+1]=ESCSr_out
            
            MCUE1[tt,:]=MCUE1_out #np.zeros((mn,t0,N));
            MCUE2[tt,:]=MCUE2_out #=np.zeros((mn,t0,N));
            ELCeq[tt]=ELCeq_out #=np.zeros((t0,N));
            ESCeq[tt]=ESCeq_out #np.zeros((t0,N));
            saturation_ratio[:,tt,:]=saturation_ratio_out #=np.zeros((mn,t0,N)); #
            gC[:,tt,:]=gC_out #np.zeros((mn,t0,N)); #growth rate at each time step
            gN[:,tt,:]=gN_out
            g[:,tt,:]=g_out #=np.zeros((mn,t0,N));
            Rm[:,tt,:]=Rm_out #=np.zeros((mn,t0,N)); #maintenance respiration at each time step
            Rg[:,tt,:]=Rg_out #=np.zeros((mn,t0,N)); #growth respiration at each time step
            Rd[tt,:]=Rd_out
            Ro[:,tt,:]=Ro_out #=np.zeros((mn,t0,N));
            RC[tt,:]=RC_out #=np.zeros((mn,t0,N));
            RL[tt,:]=RL_out #=np.zeros((mn,t0,N));
            RS[tt,:]=RS_out #=np.zeros((mn,t0,N));
            
            dGL[tt]=dGL_out #np.zeros((t0,N));
            dLMC[tt]=dLMC_out #=np.zeros((t0,N));
            dLSC[tt]=dLSC_out #=np.zeros((t0,N));
            dSAC[tt]=dSAC_out #=np.zeros((t0,N));
            dSSC[tt]=dSSC_out #=np.zeros((t0,N));
            dSPC[tt]=dSPC_out #=np.zeros((t0,N));
            BAg[:,tt,:]=BAg_out #=np.zeros((mn,t0,N));
            BAd[:,tt,:]=BAd_out
            Rma[:,tt,:]=Rma_out #=np.zeros((mn,t0,N));
            Rmd[:,tt,:]=Rmd_out #=np.zeros((mn,t0,N));
            AvailCstab[tt,:]=AvailCstab_out #=np.zeros((t0,N));
            AvailCmob[tt,:]=AvailCmob_out #=np.zeros((t0,N));
            AvailCloss[tt,:]=AvailCloss_out #=np.zeros((t0,N));
            ELCg[:,tt,:]=ELCg_out #=np.zeros((mn,t0,N));
            ESCg[:,tt,:]=ESCg_out #=np.zeros((mn,t0,N));
            
            tt=tt+1
    
        np.save(outdirtmp+'T.npy',T)
        np.save(outdirtmp+'LCin.npy',LCin)
    
        np.save(outdirtmp+'BA.npy',BA)
        np.save(outdirtmp+'BD.npy',BD)
        np.save(outdirtmp+'BAg.npy',BAg)
        np.save(outdirtmp+'BAd.npy',BAd)
        np.save(outdirtmp+'ELC.npy',ELC)
        np.save(outdirtmp+'ESC.npy',ESC)
        np.save(outdirtmp+'ELCg.npy',ELCg)
        np.save(outdirtmp+'ESCg.npy',ESCg)        
    
        np.save(outdirtmp+'GL.npy',GL)
        np.save(outdirtmp+'LMC.npy',LMC)
        np.save(outdirtmp+'LSC.npy',LSC)
        np.save(outdirtmp+'LLf.npy',LLf)
    
        np.save(outdirtmp+'SAC.npy',SAC)
        np.save(outdirtmp+'SSC.npy',SSC)
        np.save(outdirtmp+'SPC.npy',SPC)
        np.save(outdirtmp+'AvailC.npy',AvailC)
        np.save(outdirtmp+'AdsorbC.npy',AdsorbC)
    
        np.save(outdirtmp+'LMN.npy',LMN)
        np.save(outdirtmp+'LSN.npy',LSN)
    
        np.save(outdirtmp+'SAN.npy',SAN)
        np.save(outdirtmp+'SSN.npy',SSN)
        np.save(outdirtmp+'SPN.npy',SPN)
        np.save(outdirtmp+'AvailN.npy',AvailN)
        np.save(outdirtmp+'AdsorbN.npy',AdsorbN)
    
        np.save(outdirtmp+'BALr.npy',BALr)
        np.save(outdirtmp+'BASr.npy',BASr)
        np.save(outdirtmp+'BDLr.npy',BDLr)
        np.save(outdirtmp+'BDSr.npy',BDSr)
    
        np.save(outdirtmp+'ELCLr.npy',ELCLr)
        np.save(outdirtmp+'ELCSr.npy',ELCSr)
        np.save(outdirtmp+'ESCLr.npy',ESCLr)
        np.save(outdirtmp+'ESCSr.npy',ESCSr)
    
        np.save(outdirtmp+'SACLr.npy',SACLr)
        np.save(outdirtmp+'SACSr.npy',SACSr)
        np.save(outdirtmp+'SSCLr.npy',SSCLr)
        np.save(outdirtmp+'SSCSr.npy',SSCSr)
        np.save(outdirtmp+'SPCLr.npy',SPCLr)
        np.save(outdirtmp+'SPCSr.npy',SPCSr)
    
        np.save(outdirtmp+'gC.npy',gC)
        np.save(outdirtmp+'gN.npy',gN)
        np.save(outdirtmp+'g.npy',g)
        np.save(outdirtmp+'RL.npy',RL)
        np.save(outdirtmp+'RS.npy',RS)
        np.save(outdirtmp+'RC.npy',RC)
        np.save(outdirtmp+'Rm.npy',Rm)
        np.save(outdirtmp+'Rg.npy',Rg)
        np.save(outdirtmp+'Ro.npy',Ro)
        np.save(outdirtmp+'Rd.npy',Rd)
        np.save(outdirtmp+'Rma.npy',Rma)
        np.save(outdirtmp+'Rmd.npy',Rmd)
    
        np.save(outdirtmp+'AvailCLr.npy',AvailCLr)
        np.save(outdirtmp+'AvailCSr.npy',AvailCSr)
        np.save(outdirtmp+'AdsorbCLr.npy',AdsorbCLr)
        np.save(outdirtmp+'AdsorbCSr.npy',AdsorbCSr)
    
        np.save(outdirtmp+'saturation_ratio.npy',saturation_ratio)
    
        np.save(outdirtmp+'MCUE1.npy',MCUE1)
        np.save(outdirtmp+'MCUE2.npy',MCUE2)

#endtime = time.time()
#print 'elapsed time:', endtime - starttime

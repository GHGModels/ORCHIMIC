#!/usr/bin/env python
import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import os
#import copy as cp
#import time
import sys
#from scipy.optimize import minimize
#import random
#import scipy.stats as st
#import xlsxwriter

def microbe_growth(OK_check,OK_print,OK_loss,OK_N,OK_vegNuptake,N,ts,RLRS_method,K_CAE,VmaxUptake_refin,Kdin,dMFTref,Kein,Kr_ref,b,dENZ,VmaxLMC,KMLC,Ea_stab,Ea_mob,EaGL,EaKM,EaLM,EaLS,VmaxSSC,KMSC,EaSS,AdsorbCmax,KdesC,KabsC,CAE,Adj_GL,Adj_LS,Adj_SA,Adj_SP,Ke_min,Ea_uptake,Ea_main,EaSA,EaSP,KM_uptakein,ELrin,pH0,T0,H0,pHs,Ts,Hs,pH0_dec,pHs_dec,KlossC,KlossN,KlossN2,LCin,LLfin,LCNin,AvailCin,AvailNin,GL,LMC,LMN,LSC,LSN,LLf,SAC,SAN,SSC,SSN,SPC,SPN,AvailC,AvailN,AdsorbC,AdsorbN,BA,BD,MFTtype,BCN,ELC,ESC,ECN,OK_control,OK_constant_CUE,OK_stab_C,T,H,pH,clay,mn,sC,sN,AvailCLr,AvailCSr,AdsorbCLr,AdsorbCSr,SACLr,SACSr,SSCLr,SSCSr,SPCLr,SPCSr,BALr,BASr,BDLr,BDSr,ELCLr,ELCSr,ESCLr,ESCSr,LMtoSS,LStoSS,SStoSP,SAtoSP,expc,efLS,efSA,efSS,efSP,KME_FOM,KME_SOM,KME_Availin):
  zerol=1.e-5       #minimum value set in the model to avoid nan or inf
  MFT_min=1.e-5     #minimum value for microbial biomass forced in the model if set OK_control='y'
  T_ref=273.15+12.  #reference temperature, constant parameter
  T_exp=T_ref 
  max_avail=0.95  #maximum fraction of C or N in Avail pool can be used, not used anymore, Kaier et al., 2015
  #minimum and mximum CAE value from Six et al., 2006
  CAE_min=0.01
  CAE_max=0.85 #Six et al., 2006

  #senstive to temperature or not for following parameters, if yes set to 'n', otherwise set to 'y'
  OK_fixed_dMFT='y'  #default: y
  OK_fixed_KM='n'    #default: n
  OK_fixed_dENZ='y'  #default: y
  OK_fixed_Kr='n'    #default: n

  if OK_check=='y':
    #make sure LLf is not nan
    if np.sum(np.isnan(LLf)):
      print 'Error: Nan in LLf'
      if OK_print=='y':
        print 'LLf',LLf
      sys.exit()
  
  #all parameters are in unit of h, so transform to with unit of time step
  VmaxUptake_refin=VmaxUptake_refin*ts  #maximum uptake rate per active microbial biomass
  Kdin            =Kdin*ts              #To set non-linear death rate, not used in this version (Kdin=0 always)
  dMFTref         =dMFTref*ts           #death rate of active microbes at reference temperature
  Kein            =Kein*ts              #maximum enzyme production coefficient
  Kr_ref          =Kr_ref*ts            #maintenance respiration coefficient of active microbes
  dENZ            =dENZ*ts              #turn over rate of enzymes
  VmaxLMC         =VmaxLMC*ts           #maximum decomposition rate of metabolic litter
  VmaxSSC         =VmaxSSC*ts           #maximum decomposition rate of structural litter
  KabsC	       	=KabsC*ts               #
  KdesC		=KdesC*ts
  KlossC          =KlossC*ts
  KlossN	=KlossN*ts
  
  AvailCin        =AvailCin*ts
  AvailNin        =AvailNin*ts
  LCin            =LCin*ts
  LNin            =LCin/LCNin
  
  LNin[LCNin>999.]=0.
  Nplant=LNin
  
  #paramters considering the efficiency of decomposed C going to Avail pool, all set to 1 in this version, so no respiration is considered during decomposition 
  efSS[efSS<efSP]=efSP[efSS<efSP]
  efSA[efSA<efSS]=efSS[efSA<efSS]
  efLS[efLS<efSS]=efSS[efLS<efSS] #efLS is applied to both LM and LS
  
  if OK_loss=='n':
    #leaching, currently set to zero
    KlossC=0.+np.zeros((N))
    KlossN=0.+np.zeros((N))
  
  KabsC=KabsC*np.exp((0.-Ea_stab)/0.008314/T)/np.exp((0.-Ea_stab)/0.008314/T_ref)
  KdesC=KdesC*np.exp((0.-Ea_mob)/0.008314/T)/np.exp((0.-Ea_mob)/0.008314/T_ref)
  
  SAtoSS=0.15+0.68*clay-SAtoSP  #Following CENTURY model see Parton et al., 1987
  
  mn=np.size(MFTtype)
  VmaxUptake_ref=np.zeros((mn,N))
  Kd=np.zeros((mn,N))
  Ke=np.zeros((mn,N))
  dMFT_ref=np.zeros((mn,N))
  Er=np.zeros((mn,2,N))	
  
  for m in np.arange(mn):
    if MFTtype[m]==0:  #cheater, just for test
      VmaxUptake_ref[m]=VmaxUptake_refin
      Kd[m]=Kdin
      dMFT_ref[m]=dMFTref
      Ke[m]=0.
      Er[m,0]=0.
      Er[m,1]=0.
      ELC[m]=0.
      ESC[m]=0.
    elif MFTtype[m]==1: #2:
      VmaxUptake_ref[m]=VmaxUptake_refin#/0.75
      Kd[m]=Kdin#/0.75
      dMFT_ref[m]=dMFTref#/0.75
      Ke[m]=Kein
      Er[m,0]=ELrin[MFTtype[m]]
      Er[m,1]=1.-Er[m,0]
    elif MFTtype[m]==2: #5:
      VmaxUptake_ref[m]=VmaxUptake_refin#/0.75*0.5
      Kd[m]=Kdin#/0.75*0.5
      dMFT_ref[m]=dMFTref#/0.75*0.5
      Ke[m]=Kein
      Er[m,0]=ELrin[MFTtype[m]]
      Er[m,1]=1.-Er[m,0]
    if MFTtype[m]==3: #7:
      VmaxUptake_ref[m]=VmaxUptake_refin
      Kd[m]=Kdin
      dMFT_ref[m]=dMFTref
      Ke[m]=Kein
      Er[m,0]=ELrin[MFTtype[m]]
      Er[m,1]=1.-Er[m,0]
  
  
  #calculate Eeq preparing for enzymatic decompostion
  ELCeq=np.zeros((1,N))
  ESCeq=np.zeros((1,N))
  
  #ELC and ESC have a shape of (mn,3,N)
  for m in np.arange(mn):
    ELCeq=ELCeq+ELC[m]   #+ELN[m]*adjE+ELP[m]*adjE
    ESCeq=ESCeq+ESC[m]   #+ESN[m]*adjE+ESP[m]*adjE
  if OK_check=='y':
    if np.sum(ESCeq<0):
      print 'Error: negative ESCeq'
      if OK_print=='y':
        print 'ESC',ESC
      sys.exit()

  LM_fraction=0.85-0.018*LLfin*LCNin #/MMC*MMN  #LCN is molar ratio
  #To avoid lignin fraction being allocated to LM
  LM_fraction[LM_fraction>(1.-LLfin)]=(1.-LLfin)[LM_fraction>(1.-LLfin)]
  #LM_fraction[LM_fraction>1.]=1.
  LLfout=(LSC*LLf+LCin*LLfin)/(LSC+LCin*(1.-LM_fraction))
  #If LS=0, set LLf=0
  LLfout[(LSC+LCin*(1.-LM_fraction))<=0.]=0.
  LMCin=LCin*LM_fraction
  LSCin=LCin-LMCin
  
  LSNin=LNin*(1.-LM_fraction)/150.    #CN ratio of structure litter is assummed to be 150 in Parton et al., 1987
  LSNin[LSNin>LNin]=LNin[LSNin>LNin]  #When N need to keep C/N ratio of LS being 150 is larger what in input, then set all the N is allocated to LS
  LMNin=LNin-LSNin
  
  if OK_fixed_KM=='y':
    KM_adj=np.exp(-EaKM/0.008314/T_exp)/np.exp(-EaKM/0.008314/T_ref)
  else:
    KM_adj=np.exp(-EaKM/0.008314/T)/np.exp(-EaKM/0.008314/T_ref)
  
  #calculate decompostion
  #mositure function used in ORCHIDEE (Guenet et al., 2016) see ORCHDEE /modipsl/modeles/MICT/src_parameters/constantes_var.f90 moist_coeff, src_stomate/stomate_litter.f90 moistfunc_result
  control_H_tmp=-1.1*H**2+2.4*H-0.29
  control_H=np.max([0.25,np.min([1,control_H_tmp])])
  #pH function, Wang et al., 2012
  control_pH=np.exp(-(pH-pH0_dec)**2/pHs_dec**2)
  #clay function from Parton et al., 1987
  control_clay=1.-0.75*clay
  control_lignin=np.exp(-3.0*LLfout)
  
  GLeq_adj  = 1. #GL/(KMLC+GL)*Adj_GL*np.exp(-EaGL/0.008314/T)/np.exp(-EaGL/0.008314/T_ref)*control_H*control_pH*control_clay
  LMCeq_adj = 1. #LMC/(KMLC+LMC)*1.*np.exp(-EaLM/0.008314/T)/np.exp(-EaLM/0.008314/T_ref)*control_H*control_pH
  LSCeq_adj = 1. #LSC/(KMLC+LSC)*1./Adj_LS*np.exp(-EaLS/0.008314/T)/np.exp(-EaLS/0.008314/T_ref)*control_H*control_pH*control_lignin
  if OK_check=='y':
    if np.sum(np.isnan(LSCeq_adj)):
      print 'Error: nan in LSCeq_adj'	
      if OK_print=='y':
        print '1./Adj_LS',1./Adj_LS
        print 'np.exp(-EaLS/0.008314/T)/np.exp(-EaLS/0.008314/T_ref)',np.exp(-EaLS/0.008314/T)/np.exp(-EaLS/0.008314/T_ref)
        print 'control_H',control_H
        print 'control_pH',control_pH
        print 'control_lignin',control_lignin
        print 'LLfout',LLfout
      sys.exit()	
  
  GLeq  = GL*GLeq_adj
  LMCeq = LMC*LMCeq_adj
  LSCeq = LSC*LSCeq_adj
  LCeq  = GLeq+LMCeq+LSCeq
  
  SACeq_adj = 1. #SAC/(KMSC+SAC)*Adj_SA*np.exp(-EaSA/0.008314/T)/np.exp(-EaSA/0.008314/T_ref)*control_H*control_pH*control_clay
  SSCeq_adj = 1. #SSC/(KMSC+SSC)*1.*np.exp(-EaSS/0.008314/T)/np.exp(-EaSS/0.008314/T_ref)*control_H*control_pH
  SPCeq_adj = 1. #SPC/(KMSC+SPC)*1./Adj_SP*np.exp(-EaSP/0.008314/T)/np.exp(-EaSP/0.008314/T_ref)*control_H*control_pH
  
  SACeq=SAC*SACeq_adj
  SSCeq=SSC*SSCeq_adj
  SPCeq=SPC*SPCeq_adj
  SCeq=SACeq+SSCeq+SPCeq
  
  #decomposition of fresh organic matter
  #GL is not used in this version, GL represented a FOM pool that can be more easily decomposed
  dGL=GL/(KMLC*KM_adj+LCeq+ELCeq)*VmaxLMC*Adj_GL*np.exp(-EaGL/0.008314/T)/np.exp(-EaGL/0.008314/T_ref)*control_H*control_pH*control_clay*ELCeq
  tmp=GL #+GLin
  dGL[dGL>tmp]=tmp[dGL>tmp]
  
  dLMC=LMC/(KMLC*KM_adj+LCeq+ELCeq)*VmaxLMC*np.exp(-EaLM/0.008314/T)/np.exp(-EaLM/0.008314/T_ref)*control_H*control_pH*ELCeq
  tLMC=LMC+LMCin
  dLMC[dLMC>tLMC]=tLMC[dLMC>tLMC]
  
  dLMN=dLMC*(LMN+LMNin)/tLMC
  dLMN[LMC<=0.]=0.
  
  dLSC=LSC/(KMLC*KM_adj+LCeq+ELCeq)*VmaxLMC/Adj_LS*np.exp(-EaLS/0.008314/T)/np.exp(-EaLS/0.008314/T_ref)*control_H*control_pH*control_lignin*ELCeq
  tLSC=LSC+LSCin
  dLSC[dLSC>tLSC]=tLSC[dLSC>tLSC]
  
  dLSN=dLSC*(LSN+LSNin)/tLSC
  dLSN[LSC<=0.]=0.
  
  if OK_check=='y':
    if np.sum(np.isnan(LCeq)):
      print 'Error: nan in LCeq'
      print 'GLeq',GLeq
      print 'LMCeq',LMCeq
      print 'LSCeq',LSCeq
      print 'LSC',LSC
      print 'LSCeq_adj',LSCeq_adj
      sys.exit()
    if np.sum(np.isnan(dLMC)):
      print 'Error: nan in dLMC'
      print 'dLCeq',dLCeq
      print 'LMCeq',LMCeq
      print 'LCeq',LCeq
      print 'LMCeq_adj',LMCeq_adj
      sys.exit()
  
  #decomposition of soil organic carbon
  dSAC=SAC/(KMSC*KM_adj+SCeq+ESCeq)*VmaxSSC*Adj_SA*np.exp(-EaSA/0.008314/T)/np.exp(-EaSA/0.008314/T_ref)*control_H*control_pH*control_clay*ESCeq
  dSAC[dSAC>SAC]=SAC[dSAC>SAC] #only true when time step is relative small
  
  dSSC=SSC/(KMSC*KM_adj+SAC+SSC+SPC+ESCeq)*VmaxSSC*np.exp(-EaSS/0.008314/T)/np.exp(-EaSS/0.008314/T_ref)*control_H*control_pH*ESCeq
  dSSC[dSSC>SSC]=SSC[dSSC>SSC]
  
  dSPC=SPC/(KMSC*KM_adj+SAC+SSC+SPC+ESCeq)*VmaxSSC/Adj_SP*np.exp(-EaSP/0.008314/T)/np.exp(-EaSP/0.008314/T_ref)*control_H*control_pH*ESCeq
  dSPC[dSPC>SPC]=SPC[dSPC>SPC]
  
  dSAN=dSAC*SAN/SAC
  dSAN[SAC<=0.]=0.
  
  dSSN=dSSC*SSN/SSC
  dSSN[SSC<=0.]=0.
  
  dSPN=dSPC*SPN/SPC
  dSPN[SPC<=0.]=0.
  
  if OK_check=='y':
    if np.sum(dSAC<0):
      print 'Error: negative dSAC'
      print 'dSCeq',dSCeq
      print 'SACeq',SACeq
      print 'SCeq',SCeq
      print 'SACeq_adj',SACeq_adj
      sys.exit()
  
  #update AvailC pool for MFT's uptake
  CdecLM=(dLMC*(1-LMtoSS)+dGL)*efLS
  CdecLS=(dLSC*(1-LLfout)*(1.-LStoSS)+dLSC*LLfout*0.3)*efLS
  CdecSA=(dSAC*(1.-SAtoSS-SAtoSP))*efSA
  CdecSS=(dSSC*(1.-SStoSP))*efSS
  CdecSP=dSPC*efSP
  
  NdecLM=dLMN*(1-LMtoSS)
  NdecLS=dLSN*(1-LLfout)*(1.-LStoSS)+dLSN*LLfout*0.3
  NdecSA=dSAN*(1.-SAtoSS-SAtoSP)
  NdecSS=dSSN*(1.-SStoSP)
  NdecSP=dSPN
  
  AvailCdec=CdecLM+CdecLS+CdecSA+CdecSS+CdecSP
  AvailCdecLr=(CdecLM+CdecLS+CdecSA*SACLr+CdecSS*SSCLr+CdecSP*SPCLr)/AvailCdec
  AvailCdecLr[AvailCdec==0]=0.
  
  AvailNdec=(NdecLM+NdecLS+NdecSA+NdecSS+NdecSP)*(1.-KlossN2)
  
  if OK_check=='y':
    if np.sum(AvailCdec<0):
      print 'Error: negative value in AvailCdec'
      print 'CdecLM',CdecLM
      print 'CdecLS',CdecLS
      print 'CdecSA',CdecSA
      print 'CdecSS',CdecSS
      print 'CdecSP',CdecSP
      print 'dSAC',dSAC
      print 'SAtoSS',SAtoSS
      print 'SAtoSP',SAtoSP
      sys.exit()
    if np.sum(np.isnan(AvailCdec)):
      print 'Error: nan in AvailCdec'	
      print 'CdecLM',CdecLM
      print 'CdecLS',CdecLS
      print 'dGL',dGL
      print 'dLMC',dLMC
      print 'dLSC',dLSC
      print 'CdecSA',CdecSA
      print 'CdecSS',CdecSS
      print 'CdecSP',CdecSP
      sys.exit()
  
  if OK_fixed_dMFT=='y':
    dMFT_adj_T=np.exp(-Ea_uptake/0.008314/T_exp)/np.exp(-Ea_uptake/0.008314/T_ref)
  else:
    dMFT_adj_T=np.exp(-Ea_uptake/0.008314/T)/np.exp(-Ea_uptake/0.008314/T_ref)
  #dMFT_adj_Biomass=Kd*TMFT
  #dMFT_adj_Biomass[dMFT_adj_Biomass>dMFT_max]=dMFT_max
  #dMFT=dMFT_adj_Biomass*dMFT_adj_T
  #dMFT_adj_Biomass=1.+np.sum(BA+BD,axis=0)*Kd
  #dMFT=dMFT_ref*dMFT_adj_T*dMFT_adj_Biomass #*np.sum(BA,axis=0)*Kd
  dMFT=dMFT_ref*dMFT_adj_T #*np.sum(BA+BD,axis=0)
  BAd=BA*dMFT
  BAd[BAd>BA]=BA[BAd>BA]
  
  if OK_check=='y':
    if np.sum(BAd>BA):
      print 'Warnings: BAd>BA'
      print 'BAd',BAd
      print 'BA',BA
      sys.exit()
  else:
    BAd[BAd>BA]=BA[BAd>BA]
  BANd=BAd/BCN[MFTtype]
  
  AvailCmicd=0.
  AvailNmicd=0.
  AvailCmicdL=0.
  for m in np.arange(mn):
    AvailCmicd=AvailCmicd+BAd[m]*sC
    AvailNmicd=AvailNmicd+BANd[m]*sN
    AvailCmicdL=AvailCmicdL+BAd[m]*BALr[m]*sC
  
  #AvailC_tmp=AvailC+AvailCdec+AvailCmicd
  AvailC_tmp=AvailC+AvailCmicd+AvailCin
  AvailN_tmp=AvailN+AvailNdec+AvailNmicd+AvailNin
  if OK_check=='y':
    if np.sum(np.isnan(AvailC_tmp)):
      print 'Error: nan in AvailC_tmp'
      print 'AvailC',AvailC
      print 'AvailCmob',AvailCmob
      print 'AvailCdec',AvailCdec
      print 'AvailCmicd',AvailCmicd
      sys.exit()
    if np.sum(AvailC_tmp<0.):
      print 'Error: negative AvailC_tmp'
      print 'AvailC_tmp',AvailC_tmp
      print 'AvailC',AvailC
      print 'AvailCmob',AvailCmob
      print 'AvailCdec',AvailCdec
      print 'AvailCmicd',AvailCmicd
      sys.exit()
  #AvailCLr_tmp=(AvailC*AvailCLr+CdecLM+CdecLS+CdecSA*SACLr+CdecSS*SSCLr+CdecSP*SPCLr+AvailCmicdL)/AvailC_tmp
  AvailCLr_tmp=(AvailC*AvailCLr+AvailCmicdL+AvailCin)/AvailC_tmp
  AvailCLr_tmp[AvailC_tmp<=0.]=0.
  
  #update MFT pools and calculate respiration
  tmp1=-(T-T0)**2
  tmp2=Ts**2
  tmp3=tmp1/tmp2
  fT=np.exp(tmp3)
  
  fT0=np.exp(-Ea_uptake/0.008314/T)/np.exp(-Ea_uptake/0.008314/T_ref)	
  
  tmp1=-(H-H0)**2
  tmp2=Hs**2
  tmp3=tmp1/tmp2
  fH=np.exp(tmp3)
  
  tmp1=-(pH-pH0)**2
  tmp2=pHs**2
  tmp3=tmp1/tmp2
  fpH=np.exp(tmp3)
  
  VmaxUptake=VmaxUptake_ref*fT*fH*fpH*fT0
  VmaxUptakeN=VmaxUptake/BCN[MFTtype]
  fBA=BA
  
  KM_uptake_adj=np.zeros((mn,N))+VmaxUptake_ref/VmaxUptake_refin
  KMC_uptake=KM_uptakein*KM_uptake_adj*KM_adj
  KMN_uptake=KM_uptakein*KM_uptake_adj*KM_adj/BCN[MFTtype]
  
  saturation_ratio=np.zeros((mn,N))
  saturation_ratio_N=np.zeros((mn,N)) #actually not used in this version
  for i in range(mn):
    tmpC=np.zeros((N))
    tmpN=np.zeros((N))
    for j in range(mn):
      tmpC=tmpC+BA[j]/KMC_uptake[j]*KMC_uptake[i]
      tmpN=tmpN+BA[j]/KMN_uptake[j]*KMN_uptake[i]
    
    saturation_ratio[i]=AvailC/(KMC_uptake[i]+AvailC+tmpC)
    saturation_ratio_N[i]=AvailN/(KMN_uptake[i]+AvailN+tmpN)
  
  tAvailC=AvailC_tmp+AvailCdec
  tAvailN=AvailN_tmp
  Uptake=VmaxUptake*fBA*saturation_ratio
  #UptakeN=VmaxUptakeN*fBA*saturation_ratio_N
  UptakeN=Uptake*tAvailN/tAvailC
  
  if OK_check=='y':
    if np.sum(np.isnan(saturation_ratio)):
      print 'Error: nan in saturation_ratio'
      print 'AvailC_tmp',AvailC_tmp
      print 'KM_uptake',KM_uptake
      print 'KM_uptake_adj',KM_uptake_adj
      sys.exit()
    if np.sum(saturation_ratio>1.):
      print 'Error:saturation_ratio is larger than 1'
      print 'AvailC[AvailC<0.]',AvailC[AvailC<0.]
      print 'KMC_uptake[KMC_uptake<0.]',KMC_uptake[KMC_uptake<0.]
      print 'np.sum(BA,axis=0)[np.sum(BA,axis=0)<0.]',np.sum(BA,axis=0)[np.sum(BA,axis=0)<0.]
      print 'AvailN[AvailN<0.]',AvailN[AvailN<0.]
      sys.exit()
    if np.sum(np.isnan(Uptake)):
      print 'Error: nan in Uptake'
      print 'VmaxUptake',VmaxUptake
      print 'fBA',fBA
      print 'saturation_ratio',saturation_ratio
      sys.exit()
  
  TUptake=np.zeros((1,N))
  TUptakeN=np.zeros((1,N))
  
  for m in np.arange(mn):
    TUptake=TUptake+Uptake[m]
    TUptakeN=TUptakeN+UptakeN[m]
  	
  #calculate adjust uptake because total uptake should not exceed what is available
  TUptake_adj=TUptake*1.0
  TUptake_adj[TUptake_adj>tAvailC]=tAvailC[TUptake_adj>tAvailC] #uptake cannot exceed what is available
  factor_adj=TUptake_adj/TUptake
  factor_adj[TUptake==0.]=0.
  Uptake_adj=Uptake*factor_adj
  
  if RLRS_method==1:
    TUptake_adj_Lr=(AvailC_tmp*AvailCLr_tmp+AvailCdec*AvailCdecLr)/tAvailC
  elif RLRS_method==2:
    TUptake_adj_Lr=np.zeros(TUptake_adj.shape)
    TUptake_adj_Lr[TUptake_adj<=AvailCdec]=AvailCdecLr[TUptake_adj<=AvailCdec]
    TUptake_adj_Lr[TUptake_adj>AvailCdec]=((AvailCdec*AvailCdecLr+(TUptake_adj-AvailCdec)*AvailCLr_tmp)/TUptake_adj)[TUptake_adj>AvailCdec]
  else:
    print 'Error: RLRS_method has not been defined'
    sys.exit()
  TUptake_adj_Lr[TUptake_adj<=0.]=0.
  
  TUptakeN_adj=TUptakeN*1.0
  TUptakeN_adj[TUptakeN_adj>tAvailN]=tAvailN[TUptakeN_adj>tAvailN] #uptake cannot exceed what is available
  factor_adj_N=TUptakeN_adj/TUptakeN
  factor_adj_N[TUptakeN<=0.]=0.
  UptakeN_adj=UptakeN*factor_adj_N
  
  if OK_check=='y':
    if np.sum(np.isnan(Uptake_adj)):
      print 'Error: nan in Uptake_adj'
      print 'TUptake_adj',TUptake_adj
      print 'Uptake',Uptake
      print 'TUptake',TUptake
      sys.exit()	
  
  if OK_fixed_Kr=='y':
    Kr_adj=np.exp(-Ea_main/0.008314/T_exp)/np.exp(-Ea_main/0.008314/T_ref)
  else:
    Kr_adj=np.exp(-Ea_main/0.008314/T)/np.exp(-Ea_main/0.008314/T_ref)
  Kr=Kr_ref*Kr_adj		
  
  Rma=BA*Kr
  tmp=BA+Uptake_adj-BAd
  Rma[Rma>tmp]=tmp[Rma>tmp]
  
  if OK_check=='y':
    if np.sum(np.isnan(TUptake_adj)):
      print 'Error: nan in TUptake_adj'
      print 'TUptake',TUptake
      print 'Uptake',Uptake
      sys.exit()
    if np.sum(Rma<0):
      print 'Error: negative Rma'
      print 'BA',BA
      print 'Kr',Kr
      sys.exit()

  if OK_constant_CUE=='y':
    CAEtmp=CAE*1.
  elif OK_constant_CUE=='n':
    CAEtmp=CAE-K_CAE*(T-T_ref)
    CAEtmp[CAEtmp<CAE_min]=CAE_min
    CAEtmp[CAEtmp>CAE_max]=CAE_max
  
  Cgrowth=Uptake_adj-Rma
  Cgrowth[Cgrowth<=0.]=0.
  gC=Cgrowth/BA*CAEtmp
  gC[BA<=0.]=0.
  
  Ngrowth=UptakeN_adj
  gN=Ngrowth/(BA/BCN[MFTtype])
  
  if OK_N=='y':
    g=np.where(gC>gN,gN,gC)
  elif OK_N=='n':
    g=gC
  
  BAg=BA*g
  Rg=BAg*(1-CAEtmp)/CAEtmp	
  
  if OK_N=='y':
    AvailNuptake_rlease=np.sum(Ngrowth-BAg/BCN[MFTtype],axis=0)
  elif OK_N=='n':
    AvailNuptake_rlease=np.zeros((N))		
  
  Ro=Cgrowth-BAg-Rg
  Ro[Ro<=0.]=0.
  respLM=(dLMC*(1-LMtoSS)+dGL)*(1-efLS)
  respLS=(dLSC*(1-LLfout)*(1.-LStoSS)+dLSC*LLfout*0.3)*(1-efLS)
  respSA=(dSAC*(1.-SAtoSS-SAtoSP))*(1-efSA)
  respSS=(dSSC*(1.-SStoSP))*(1-efSS)
  respSP=dSPC*(1-efSP)
  #Rd is respiration from decomposition, but actually not considered in this version, so Rd=0.
  Rd=respLM+respLS+respSA+respSS+respSP
  RdL=respLM+respLS+respSA*SACLr+respSS*SSCLr+respSP*SPCLr
  
  BAmain=Rma-Uptake_adj
  BAmain[BAmain<0.]=0.
  tmp=BA-BAd
  BAmain[BAmain>tmp]=tmp[BAmain>tmp]
  
  BAmainN=BAmain/BCN[MFTtype]
  
  BAtoD=(1-saturation_ratio)*Kr*BA
  BDtoA=saturation_ratio*Kr*BD
  
  tmp=BA-BAd-BAmain
  BAtoD[BAtoD>tmp]=tmp[BAtoD>tmp]
  
  if OK_check=='y':
    if np.sum(BAtoD<-zerol):
      print 'Warning: negative BAtoD'
      if np.sum(tmp<-1e-5):
        print 'BA[tmp<-1e-5]',BA[tmp<-1e-5]
        print 'BAg[tmp<-1e-5]',BAg[tmp-1e-5]
        print 'BAd[tmp<-1e-5]',BAd[tmp<-1e-5]
        print 'Rma[tmp<-1e-5]',Rma[tmp<-1e-5]
      print 'BAtoD[BAtoD<-1e-5]',BAtoD[BAtoD<-1e-5]
      print 'tmp[[BAtoD<-1e-5]]',tmp[[BAtoD<-1e-5]]
      print 'Kra[Kr<0.]',Kr[Kr<0.]
      print 'BA[BA<0.]',BA[BA<0.]
      print 'saturation_ratio[saturation_ratio<0.]',saturation_ratio[saturation_ratio<0.]
      print 'saturation_ratio[saturation_ratio>1.]',saturation_ratio[saturation_ratio>1.]
      sys.exit()
    else:
      BAtoD[BAtoD<0.]=0.
  else:
  	BAtoD[BAtoD<0.]=0.
  
  Rmd=BD*b*Kr
  tmp=BD
  Rmd[Rmd>tmp]=tmp[Rmd>tmp]
  
  if OK_check=='y':
    if np.sum(Rmd<-zerol):
      print 'Error: negative Rmd'
      print 'BD[Rmd<-zerol]',BD[Rmd<-zerol]
      print 'np.sum(BD<0.)',np.sum(BD<0.)
      print 'b',b
      print 'Kr',Kr
      print 'tmp[tmp<0.]',tmp[tmp<0.]
      sys.exit()
    else:
      Rmd[Rmd<0.]=0.
  else:
    Rmd[Rmd<0.]=0.
  
  tmp=BD-Rmd
  BDtoA[BDtoA>tmp]=tmp[BDtoA>tmp]
  
  BDmain=Rmd
  Rm=Rma+Rmd
  
  BDmainN=BDmain/BCN[MFTtype]
  
  Tresp=np.sum(Rm+Rg+Ro,axis=0)+Rd
 
  #MCUE1 is the CUE used in the paper
  MCUE1=np.sum(BAg,axis=0)/(Tresp+np.sum(BAg,axis=0))
  MCUE1[(Tresp+np.sum(BAg,axis=0))==0]=0.
  MCUE2=np.sum(BAg,axis=0)/TUptake_adj
  MCUE2[TUptake_adj==0.]=0.
  
  RC=Tresp #Rg+Rm+Ro
  CAmain=BAmain
  NAmain=BAmainN
  
  RL=np.sum((Rg+Rma-CAmain+Ro)*TUptake_adj_Lr+Rmd*BDLr+CAmain*BALr,axis=0)+RdL
  RS=RC-RL
  
  #update E pools
  ELC_factor1=LCeq/(KME_FOM*KM_adj+LCeq+ELCeq)		
  ESC_factor1=SCeq/(KME_SOM*KM_adj+SCeq+ESCeq)
  
  KME_Avail_adj=KM_uptake_adj
  KME_Avail=KME_Availin*KME_Avail_adj*KM_adj
  EC_factor2=np.zeros((mn,N))
  for i in range(mn):
    tmpC=np.zeros((N))
    #tmpN=np.zeros((N))
    for j in range(mn):
      tmpC=tmpC+BA[j]/KME_Avail[j]*KME_Avail[i]
      #tmpN=tmpN+BA[j]/KMN_uptake[j]*KMN_uptake[i]
    
    EC_factor2[i]=1.-AvailC/(KME_Avail[i]+AvailC+tmpC)
  
  
  Ke_min_2D=np.zeros((mn,N))
  Ke_min_2D[:,:]=Ke_min
  
  ELC_factor=(ELC_factor1*EC_factor2)**expc #expc=0.5. Sinsabaugh et al., 2013
  ESC_factor=(ESC_factor1*EC_factor2)**expc 
  
  ELC_factor[ELC_factor<Ke_min_2D]=Ke_min_2D[ELC_factor<Ke_min_2D]
  ESC_factor[ESC_factor<Ke_min_2D]=Ke_min_2D[ESC_factor<Ke_min_2D]
  
  KeCL=Ke*Er[:,0]*ELC_factor
  KeCS=Ke*Er[:,1]*ESC_factor
  
  BAenz=BA*(KeCL+KeCS)		
  
  BAenz_ori=BAenz*1.
  tmp=BA-BAd-BAmain-BAtoD
  BAenz[BAenz>tmp]=tmp[BAenz>tmp]
  KeCL=KeCL*BAenz/BAenz_ori
  KeCS=KeCS*BAenz/BAenz_ori
  
  KeCL[BAenz<=0.]=0.
  KeCS[BAenz<=0.]=0.
  
  ELCg=KeCL*BA
  ESCg=KeCS*BA
  
  if OK_fixed_dENZ=='y':
    dENZ_adj=np.exp(-Ea_uptake/0.008314/T_exp)/np.exp(-Ea_uptake/0.008314/T_ref)
  else:
    dENZ_adj=np.exp(-Ea_uptake/0.008314/T)/np.exp(-Ea_uptake/0.008314/T_ref)
  ELCd=ELC*dENZ*dENZ_adj
  ESCd=ESC*dENZ*dENZ_adj
  
  tmp=ELC #+ELCg
  ELCd[ELCd>tmp]=tmp[ELCd>tmp]
  tmp=ESC #+ESCg
  ESCd[ESCd>tmp]=tmp[ESCd>tmp]
  
  ELCout=ELC+ELCg-ELCd
  #ELNout=ELN+ELNg-ELNd
  #ELPout=ELP+ELPg-ELPd
  ESCout=ESC+ESCg-ESCd
  #ESNout=ESN+ESNg-ESNd
  #ESPout=ESP+ESPg-ESPd
  
  if OK_check=='y':
    if np.sum(ELCout<0.) or np.sum(ESCout<0.):
      print 'Warning: ELCout or ESCout is negative'
      print 'ELCout',ELCout
      print 'ESCout',ELSout
      sys.exit()
  else:
    ELCout[ELCout<0.]=0.
    ESCout[ESCout<0.]=0.
  
  ELCLrout=(ELC*ELCLr+ELCg*BALr-ELCd*ELCLr)/ELCout
  #ELNLrout=(ELN*ELNLr+ELNg*BALr-ELNd*ELNLr)/ELNout
  #ELPLrout=(ELP*ELPLr+ELPg*BALr-ELPd*ELPLr)/ELPout
  ESCLrout=(ESC*ESCLr+ESCg*BALr-ESCd*ESCLr)/ESCout
  #ESNLrout=(ESN*ESNLr+ESNg*BALr-ESNd*ESNLr)/ESNout
  #ESPLrout=(ESP*ESPLr+ESPg*BALr-ESPd*ESPLr)/ESPout
  ELCLrout[ELCout<=0.]=0.
  ESCLrout[ESCout<=0.]=0.
  ELCSrout=1.-ELCLrout
  #ELNSrout=1.-ELNLrout
  #ELPSrout=1.-ELPLrout
  ESCSrout=1.-ESCLrout
  #ESNSrout=1.-ESNLrout
  #ESPSrout=1.-ESPLrout
  
  BAout=BA+BAg-BAmain-BAenz-BAd-BAtoD+BDtoA
  BDout=BD-BDmain+BAtoD-BDtoA		
  
  BALrout=(BA*BALr+BAg*TUptake_adj_Lr-BAmain*BALr-BAenz*BALr-BAd*BALr-BAtoD*BALr+BDtoA*BDLr)/BAout
  #BASrout=(BA*BASr+BAg*AvailCSr-BAmain*BASr-BAenz*BASr-BAd*BASr-BAtoD*BASr+BDtoA*BDSr)/BAout
  BALrout[BAout<=0.]=0.
  BALrout[BALrout<0.]=0.
  BALrout[BALrout>1.]=1.
  BASrout=1.- BALrout
  
  BDLrout=(BD*BDLr-BDmain*BDLr+BAtoD*BALr-BDtoA*BDLr)/BDout
  #BDSrout=(BD*BDSr-BDmain*BDSr+BAtoD*BASr-BDtoA*BDSr)/BDout
  BDLrout[BDout<=0.]=0.
  BDLrout[BDLrout<0.]=0.
  BDLrout[BDLrout>1.]=1.
  BDSrout=1.- BDLrout
  
  if OK_check=='y':
    if np.sum(BAout<-zerol) or np.sum(BDout<-zerol):
      print 'Warning: BAout or BDout is negative'
      print 'BAout[BAout<0.]',BAout[BAout<0.]
      print 'BDout[BDout<0.]',BDout[BDout<0.]
      sys.exit()
    else:
      BAout[BAout<0.]=0.
      BDout[BDout<0.]=0.
    if np.sum(BALrout<0.):
      print 'Error: Negative BALrout'
      print 'BAout[BALrout<0.]',BAout[BALrout<0.]
      print '(BA*BALr+BAg*AvailCLr-BAmain*BALr-BAenz*BALr-BAd*BALr-BAtoD*BALr+BDtoA*BDLr)[BALrout<0.]',(BA*BALr+BAg*AvailCLr-BAmain*BALr-BAenz*BALr-BAd*BALr-BAtoD*BALr+BDtoA*BDLr)[BALrout<0.]
      sys.exit()
  else:
    BAout[BAout<0.]=0.
    BDout[BDout<0.]=0.
  
  if OK_check=='y':
    if np.sum(np.isnan(BAout)):
      print 'Error: NAN in BAout'
      print 'BAout',BAout
      print 'BA',BA
      print 'BAg',BAg
      print 'BAmain',BAmain
      print 'BAenz',BAenz
      print 'BAd',BAd
      print 'BAtoD',BAtoD
      print 'BDtoA',BDtoA
      sys.exit()
  else:
      BAout[BAout<0.]=0.
      BDout[BDout<0.]=0.
  
  #update Avail pool
  AvailC2_tmp=AvailC_tmp+AvailCdec-TUptake_adj
  AvailC2Lr_tmp=(AvailC_tmp*AvailCLr_tmp+AvailCdec*AvailCdecLr-TUptake_adj*TUptake_adj_Lr)/AvailC2_tmp
  AvailC2Lr_tmp[AvailC2_tmp==0.]=0
  if np.sum(np.isnan(AvailC2_tmp)):
    print 'Error in lin 993 for AvailC2_tmp'
    print 'AvailC_tmp',AvailC_tmp
    print 'AvailCLr_tmp',AvailCLr_tmp
    print 'AvailCdec',AvailCdec
    print 'AvailCdecLr',AvailCdecLr
    print 'TUptake_adj',TUptake_adj
    print 'TUptake_adj_Lr',TUptake_adj_Lr
    print 'AvailC2_tmp',AvailC2_tmp
    sys.exit()
  
  AvailCmob=KdesC*AdsorbC/AdsorbCmax
  AvailCmob[AvailCmob>AdsorbC]=AdsorbC[AvailCmob>AdsorbC]
  
  AvailNmob=AvailCmob*AdsorbN/AdsorbC
  AvailNmob[AdsorbC==0.]=0.
  #print AvailNmob,AvailCmob,AdsorbN,AdsorbC
  
  AvailCstab=AvailC*KabsC*(1-AdsorbC/AdsorbCmax)
  AvailCstab[AvailCstab<0.]=0.
  tmp=AvailC2_tmp #-TUptake_adj
  AvailCstab[AvailCstab>tmp]=tmp[AvailCstab>tmp]
  
  AvailNstab=AvailCstab*AvailN/AvailC
  tmp=AvailN_tmp-TUptakeN_adj
  AvailNstab[AvailNstab>tmp]=tmp[AvailNstab>tmp]
  
  AvailCloss=AvailC*KlossC
  tmp=AvailC2_tmp-AvailCstab
  AvailCloss[AvailCloss>tmp]=tmp[AvailCloss>tmp]
  
  #AvailNloss=AvailN*KlossN
  AvailNloss2=AvailCloss*AvailN/AvailC #leach
  AvailNloss=(NdecLM+NdecLS+NdecSA+NdecSS+NdecSP)*KlossN2+AvailNloss2  #Same as Parton et al., 1987, %5 of totoal mineralization flux as volatilization and leach was negelected
  
  AvailCout=AvailC2_tmp-AvailCloss-AvailCstab+AvailCmob
  AvailCLrout=(AvailC2_tmp*AvailC2Lr_tmp-AvailCloss*AvailC2Lr_tmp-AvailCstab*AvailC2Lr_tmp+AvailCmob*AdsorbCLr)/AvailCout
  AvailCLrout[AvailCout<=0.]=0.
  AvailCSrout=1.-AvailCLrout
  
  AvailNmain=np.sum(BAmainN+BDmainN,axis=0)
  #print 'AvailNmain.shape',AvailNmain.shape
  
  AvailNout=AvailN_tmp-TUptakeN_adj-AvailNstab+AvailNmob+AvailNmain+AvailNuptake_rlease-AvailNloss2 #Here is different from that for C, because AvailN_tmp has alreay considered AvailNloss
  
  if OK_vegNuptake=='y':
    AvailNout[AvailNout>Nplant]=(AvailNout-Nplant)[AvailNout>Nplant]
  
  #print AvailN_tmp,TUptakeN_adj,AvailNstab,AvailNmob,AvailNmain,AvailNuptake_rlease,AvailNloss2
  #sys.exit() 
  
  #print 'AvailNout.shape',AvailNout.shape
  if OK_check=='y':
    if np.sum(np.isnan(AvailCLrout)):
      print 'Nan in AvailCLrout'
      sys.exit()
    if np.sum(AvailCout<0.):
      print 'Warning: AvailCout is negative'
      print 'AvailCout',AvailCout
      sys.exit()
    if np.sum(AvailNout<0.):
      print 'Warning: AvailNout is negative'
      sys.exit()
  else:
    AvailCout[AvailCout<=0.]=0.
  
  #update three Litter C pools
  GLout=GL-dGL #+GLin
  #LM_fraction=0.85-0.018*LLfin*LCN #/MMC*MMN  #LCN is molar ratio
  #if LM_fraction>1.-LLfin: #to avoid lignin being allocated to LM
  #	LM_fraction=1.-LLfin
  LMCout=LMC-dLMC+LMCin
  LSCout=LSC-dLSC+LSCin
  LMCout[LMCout<0.]=0.
  LSCout[LSCout<0.]=0.
  #LLfout=((LSC-dLSC)*LLf+LCin*LLfin)/LSCout  #all lignin go into LS
  #LLfout[LSCout<=0.]=0.	
  
  LMNout=LMN-dLMN+LMNin
  LSNout=LSN-dLSN+LSNin
  
  #update SAC pool
  SACmicd=np.zeros((1,N))
  SANmicd=np.zeros((1,N))
  SACmicdL=np.zeros((1,N))
  for m in np.arange(mn):
    SACmicd=SACmicd+BAd[m]*(1-sC)
    SANmicd=SANmicd+BANd[m]*(1-sC)
    SACmicdL=SACmicdL+BAd[m]*BALr[m]*(1-sC)
  
  SACenzd=np.zeros((1,N))
  SANenzd=np.zeros((1,N))
  SACenzdL=np.zeros((1,N))
  for m in np.arange(mn):
    SACenzd=SACenzd+(ELCd[m]+ESCd[m]) #+ELNd[m]+ESNd[m]+ELPd[m]+ESPd[m])
    SANenzd=SANenzd+(ELCd[m]+ESCd[m])/ECN[MFTtype[m]]
    SACenzdL=SACenzdL+(ELCd[m]*ELCLr[m]+ESCd[m]*ESCLr[m]) #+ELNd[m]*ELNLr[m]+ESNd[m]*ESNLr[m]+ELPd[m]*ELPLr[m]+ESPd[m]*ESPLr[m])
  
  SACout=SAC+SACmicd+SACenzd-dSAC
  SACout[SACout<0.]=0.
  SANout=SAN+SANmicd+SANenzd-dSAN
  SACLrout=(SAC*SACLr+SACmicdL+SACenzdL-dSAC*SACLr)/SACout
  SACLrout[SACout<=0.]=0.
  SACSrout=1.-SACLrout
  
  if OK_check=='y':
    if np.sum(SACout<-zerol):
      print 'Error: negative SACout'
      print 'SACout[SACout<-zerol]',SACout[SACout<-zerol]
      sys.exit()
    else:
      SACout[SACout<0.]=0.
  else:
    SACout[SACout<0.]=0.
  
  #update SSC pool
  SSCout=SSC+dLMC*LMtoSS+(dLSC*(1-LLfout)*LStoSS+dLSC*LLfout*(1-0.3))+dSAC*SAtoSS-dSSC
  SSCout[SSCout<0.]=0
  SSNout=SSN+dLMN*LMtoSS+(dLSN*(1-LLfout)*LStoSS+dLSN*LLfout*(1-0.3))+dSAN*SAtoSS-dSSN
  
  SSCLrout=(SSC*SSCLr+dLMC*LMtoSS+(dLSC*(1-LLfout)*LStoSS+dLSC*LLfout*(1-0.3))+dSAC*SAtoSS*SACLr-dSSC*SSCLr)/SSCout
  SSCLrout[SSCout<=0.]=0.
  SSCSrout=1.-SSCLrout
  
  if OK_check=='y':
    if np.sum(SSCout<-zerol):
      print 'Error: negative SSCout'
      print 'SSCout[SSCout<-zerol]',SSCout[SSCout<-zerol]
      sys.exit()
    else:
      SSCout[SSCout<0.]=0.
  else:   
    SSCout[SSCout<0.]=0.
  
  #update SPC pool
  SPCout=SPC+dSAC*SAtoSP+dSSC*SStoSP-dSPC
  SPNout=SPN+dSAN*SAtoSP+dSSN*SStoSP-dSPN		
  
  SPCLrout=(SPC*SPCLr+dSAC*SAtoSP*SACLr+dSSC*SStoSP*SSCLr-dSPC*SPCLr)/SPCout
  SPCLrout[SPCout<=0.]=0.
  SPCSrout=1.-SPCLrout
  
  if OK_check=='y':
    if np.sum(SPCout<-zerol):
      print 'Error: negative SPCout'
      print 'SPCout[SPCout<-zerol]',SPCout[SPCout<-zerol]
      sys.exit()
    else:
      SPCout[SPCout<0.]=0.
  else:   
    SPCout[SPCout<0.]=0.	
  
  AdsorbCout=AdsorbC+AvailCstab-AvailCmob
  AdsorbNout=AdsorbN+AvailNstab-AvailNmob
  AdsorbNout[AdsorbNout<0.]=0.
  
  AdsorbCLrout=(AdsorbC*AdsorbCLr+AvailCstab*AvailC2Lr_tmp-AvailCmob*AdsorbCLr)/AdsorbCout
  AdsorbCLrout[AdsorbCout<=0.]=0.
  AdsorbCSrout=1.-AdsorbCLrout #(AdsorbC*AdsorbCSr+AdsorbCabs*AvailCSr-AdsorbCdes*AdsorbCSr)/AdsorbCout
  
  if OK_check=='y':
    if np.sum(AdsorbCout<-zerol):
      print 'Error: negative AdsorbCout'
      print 'AdsorbCout[AdsorbCout<-zerol]',AdsorbCout[AdsorbCout<-zerol]
      sys.exit()
    else:
      AdsorbCout[AdsorbCout<0.]=0.
  else:   
    AdsorbCout[AdsorbCout<0.]=0.
  
  if OK_check=='y':
    C1=AvailCin+LCin+GL+LMC+LSC+SAC+SSC+SPC+AvailC+AdsorbC+np.sum((BA+BD+ELC+ESC),axis=0)
    C2=GLout+LMCout+LSCout+SACout+SSCout+SPCout+AvailCout+AdsorbCout+np.sum((BAout+BDout+ELCout+ESCout),axis=0)+np.sum(Rma+Rmd+Rg+Ro,axis=0)+AvailCloss
    print 'C=',C1,C2
    if np.sum(np.abs(C1-C2)>1e-5):
      print 'Error: C mass is not close'
      print 'np.sum(np.abs(C1-C2)>1e-5)',np.sum(np.abs(C1-C2)>1e-5)
      print '(C1-C2)[np.abs(C1-C2)>1e-5]',(C1-C2)[np.abs(C1-C2)>1e-5]
      print 'C1[np.abs(C1-C2)>1e-5]',C1[np.abs(C1-C2)>1e-5]
      print 'C2[np.abs(C1-C2)>1e-5]',C2[np.abs(C1-C2)>1e-5]
      sys.exit()

    if OK_N=='y':
      N1=LNin+LMN+LSN+SAN+SSN+SPN+AvailN+AdsorbN+np.sum((BA+BD)/BCN[MFTtype],axis=0)+np.sum((ELC+ESC)/ECN[MFTtype],axis=0)
      N2=LMNout+LSNout+SANout+SSNout+SPNout+AvailNout+AdsorbNout+np.sum((BAout+BDout)/BCN[MFTtype],axis=0)+np.sum((ELCout+ESCout)/ECN[MFTtype],axis=0)+AvailNloss
      print 'N=',N1,N2
      if np.sum(np.abs(N1-N2)>1e-6):
        print 'Error: N mass is not close'
        print 'np.sum(np.abs(N1-N22)>1e-6)',np.sum(np.abs(N1-N2)>1e-6)
        print '(N1-N2)[np.abs(N1-N2)>1e-6]',(N1-N2)[np.abs(N1-N2)>1e-6]
        print '(C1-C2)[np.abs(N1-N2)>1e-6]',(C1-C2)[np.abs(N1-N2)>1e-6]
        print 'N1[np.abs(N1-N2)>1e-6]',N1[np.abs(N1-N2)>1e-6]
        print 'N2[np.abs(N1-N2)>1e-6]',N2[np.abs(N1-N2)>1e-6]
        sys.exit()
  
  
  #	##"Check P balance"
  #	#P1=LitterPin+LitterSP+LitterMP+ActiveP+SlowP+PassiveP+AvailP+AdsorbP+np.sum((BA+BD+ELC+ESC+ELN+ESN+ELP+ESP)*rP/rC,axis=0)+AvailPin
  #	#P2=LitterSPout+LitterMPout+ActivePout+SlowPout+PassivePout+AvailPout+AdsorbPout+VEGPg+np.sum((BAout+BDout+ELCout+ESCout+ELNout+ESNout+ELPout+ESPout)*rP/rC,axis=0)+AvailPloss
  #	#if np.sum(np.abs(P2-P1)>zerol):
  #	#	print 'Error: The mass of P is not close:',P1,P2,P2-P1
  #	#	sys.exit()
  #	
    CF1=AvailCin+LCin+GL+LMC+LSC+SAC*SACLr+SSC*SSCLr+SPC*SPCLr+AvailC*AvailCLr+AdsorbC*AdsorbCLr+np.sum((BA*BALr+BD*BDLr+ELC*ELCLr+ESC*ESCLr),axis=0)
    CF2=GLout+LMCout+LSCout+SACout*SACLrout+SSCout*SSCLrout+SPCout*SPCLrout+AvailCout*AvailCLrout+AdsorbCout*AdsorbCLrout+np.sum((BAout*BALrout+BDout*BDLrout+ELCout*ELCLrout+ESCout*ESCLrout),axis=0)+RL+AvailCloss*AvailC2Lr_tmp
    if np.sum(abs(CF1-CF2)>1e-6):
      print 'Error: CF is not closed'
      print 'CF1',CF1
      print 'CF2',CF2
      print 'C',C1,C2
      print AvailCin
      print LCin
      print GL
      print LMC
      print LSC
      print SAC*SACLr
      print SSC*SSCLr
      print SPC*SPCLr
      print AvailC*AvailCLr
      print AdsorbC*AdsorbCLr
      print np.sum((BA*BALr+BD*BDLr+ELC*ELCLr+ESC*ESCLr),axis=0)
      sys.exit()
  
  if OK_control=='y':
    for m in np.arange(mn):
      if BAout[m]<MFT_min:
        BAout[m]=MFT_min
      if BDout[m]<MFT_min:
        BDout[m]=MFT_min
  
  return GLout,LMCout,LMNout,LSCout,LSNout,LLfout,SACout,SANout,SSCout,SSNout,SPCout,SPNout,AvailCout,AvailNout,BAout,BDout,ELCout,ESCout,ELCeq,ESCeq,AdsorbCout,AdsorbNout,MCUE1,MCUE2,saturation_ratio,Rd,Ro,Rg,Rm,RC,gC,gN,g,AvailCLrout,AvailCSrout,AdsorbCLrout,AdsorbCSrout,SACLrout,SACSrout,SSCLrout,SSCSrout,SPCLrout,SPCSrout,BALrout,BASrout,BDLrout,BDSrout,ELCLrout,ELCSrout,ESCLrout,ESCSrout,dGL,dLMC,dLMN,dLSC,dLSN,dSAC,dSAN,dSSC,dSSN,dSPC,dSPN,BAg,BAd,Rma,Rmd,RL,RS,AvailCstab,AvailNstab,AvailCmob,AvailNmob,AvailCloss,AvailNloss,ELCg,ESCg

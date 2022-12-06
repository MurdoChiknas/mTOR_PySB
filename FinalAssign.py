#usr/bin/env python3
from pysb import *
from pysb.macros import bind
from pysb.macros import *
from pysb.integrate import odesolve
import numpy as np
import sys
from pysb.simulator import ScipyOdeSimulator
import pylab as pl              

###Set Starting Rapamycin & Insulin Concentration
rapa_0 = sys.argv[1]
insulin_0 = sys.argv[2]


### Insulin-IRa Kf & Kd
    #https://doi.org/10.3389/fendo.2015.00107, Ligand-binding affinity at the insulin receptor isoform-A and subsequent IR-A tyrosine phosphorylation kinetics are important determinants of mitogenic biological outcomes, Harinda Rajapaksha and imageBriony E. Forbes
ins_kf = 0.636e-9
ins_kr = 1.57e-9

### Rapamycin-Related Kf & Kd 
    #https://pubs.acs.org/doi/10.1021/ja043277y, Characterization of the FKBP·Rapamycin·FRB Ternary Complex, Laura A. Banaszynski, Corey W. Liu, and Thomas J. Wandless
#RAPA-FKBP Binding
rapa_kf = 5e-9 
rapa_kr = 0.2e-9
#RAPA-mTOR/FRB Binding
rapa_mt_kr = 26e-6 
rapa_mt_kf = 0.038e-6
#RAPA/FKBP-mTOR/FRB Complex Binding
rapa_comp_kr = 12e-9 
rapa_comp_kf = 0.0833e-9

### mTOR Catalytic Activity:
    #https://doi.org/10.1038/nature25023, Mechanisms of mTORC1 activation by RHEB and inhibition by PRAS40, Yang, Haijuan; Jiang, Xiaolu; Li, Buren; Yang, Hyo J.; Miller, Meredith; Yang, Angela; Dhar, Ankita; Pavletich, Nikola P. 
#Phosphorylation of S6K by mTORC1
kcat_ps6 = 0.66
#Phosphorylation of 4EBP by mTORC1
kcat_pebp4 = 0.003

#Initiate model
Model()

#Declare Monomers
Monomer('RAPA', ['s']) #Rapamycin
Monomer('INS', ['s']) #Insulin
Monomer('FKBP', ['sRAPA', 'sTOR'])
Monomer('mTOR', ['sRHEB', 'sATP', 'kFRB'])
Monomer('S6K', ['s'])
Monomer('pS6K', ['p'])
Monomer('EBP4', ['s'])
Monomer('pEBP4', ['p'])
Monomer('ADP', ['s'])
Monomer('ATP', ['s'])
Monomer('AKT', ['s'])
Monomer('pAKT', ['p'])
Monomer('RHEB', ['kAKT'])
Monomer('pRHEB', ['p'])
Monomer('IRS', ['sINS', 'sPI3K', 'sATP', 'kAKT'])
Monomer('PI3K', ['sS6K', 'sIRS'])
Monomer('iPI3K', ['i']) #Inactivated iPI3K

#Reaction Rates
Parameter('kf_bind', 1e-5)
Parameter('kr_bind', 1e-1)
Parameter('kcat', 1e-1)
Parameter('kcat_s6', kcat_ps6)
Parameter('kcat_ebp4', kcat_pebp4)
Parameter('rapa_kfwd', rapa_kf)
Parameter('rapa_krev', rapa_kr)
Parameter('rapa_mt_kfwd', rapa_mt_kf)
Parameter('rapa_mt_krev', rapa_mt_kr)
Parameter('rapa_comp_kfwd', rapa_comp_kf)
Parameter('rapa_comp_krev', rapa_comp_kr)
Parameter('ins_kfwd', ins_kf)
Parameter('ins_krev', ins_kr)

#Initial Concentrations

    #Exogenous Molecules
Initial(RAPA(s=None), Parameter('RAPA_0', rapa_0)) #Tunable mTOR inhibition
Initial(INS(s=None), Parameter('INS_0',insulin_0)) #Tuneable mTOR activation

    #Endogenous Molecules
Initial(FKBP(sRAPA=None, sTOR=None), Parameter('FKBP_0', 2000))
Initial(mTOR(sRHEB=None, sATP=None, kFRB=None), Parameter('mTOR_0', 10000))
Initial(S6K(s=None), Parameter('S6K_0', 6000))
Initial(EBP4(s=None), Parameter('EBP4_0',5000))
Initial(ATP(s=None), Parameter('ATP_0',10000))  
Initial(AKT(s=None), Parameter('AKT_0',5000))
Initial(RHEB(kAKT=None), Parameter('RHEB_0',1000))
Initial(IRS(sINS=None, sPI3K=None, sATP=None, kAKT=None), Parameter('IRS_0',10000))
Initial(PI3K(sS6K=None, sIRS=None), Parameter('PI3K_0',5000))


#PATHWAYS:

#Insulin Binding & IRS Activation w/ ATP --> Production of phospho-AKT
bind(INS, 's', IRS, 'sINS',[ins_kfwd, ins_krev])
Rule('Active_IRS_ATP_Bind', INS(s=1)%IRS(sINS=1,sPI3K=None, sATP=None, kAKT=None)+ ATP(s=None)| INS(s=1)%IRS(sINS=1,sPI3K=None, sATP=2, kAKT=None)%ATP(s=2), *[kf_bind, kr_bind])
Rule('Active_IRS_PI3K_Bind', INS(s=1)%IRS(sINS=1,sPI3K=None, sATP=None, kAKT=None)+ PI3K(sS6K=None, sIRS=None)| INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=None, kAKT=None)%PI3K(sS6K=None, sIRS=3), *[kf_bind, kr_bind])
Rule('Active_IRS_PI3K_Bind_ATP_Bind', INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=None, kAKT=None)%PI3K(sS6K=None, sIRS=3) + ATP(s=None)| INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=2, kAKT=None)%PI3K(sS6K=None, sIRS=3)%ATP(s=2), *[kf_bind, kr_bind])
Rule('Active_IRS_ATP_Bind_PI3K_Bind', INS(s=1)%IRS(sINS=1,sPI3K=None, sATP=2, kAKT=None)%ATP(s=2)+ PI3K(sS6K=None, sIRS=None)| INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=2, kAKT=None)%ATP(s=2)%PI3K(sS6K=None, sIRS=3), *[kf_bind, kr_bind])
Rule('Activate_AKT', INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=2, kAKT=None)%ATP(s=2)%PI3K(sS6K=None, sIRS=3)+AKT(s=None) | INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=2, kAKT=4)%ATP(s=2)%PI3K(sS6K=None, sIRS=3)%AKT(s=4), *[kf_bind, kr_bind])
Rule('Phosphorylate_AKT', INS(s=1)%IRS(sINS=1,sPI3K=3, sATP=2, kAKT=4)%ATP(s=2)%PI3K(sS6K=None, sIRS=3)%AKT(s=4) >> INS(s=None) + IRS(sINS=None,sPI3K=None, sATP=None, kAKT=None) + ADP(s=None) + PI3K(sS6K=None, sIRS=None) + pAKT(p=None), *[kcat])

#PhosphoAKT-RHEB Binding --> Prodcution of Activated RHEB
bind(pAKT, 'p', RHEB, 'kAKT',[kf_bind, kr_bind])
Rule('RHEB_Activation', pAKT(p=1)%RHEB(kAKT=1) >> AKT(s=None) + pRHEB(p=None), *[kcat])

#mTOR-RAPA-FKBP Inhibitory Complex
Rule('FKBP_Inhibitory_Complex', RAPA(s=None)+ FKBP(sRAPA=None, sTOR=None) | RAPA(s=2)%FKBP(sRAPA=2, sTOR=None), *[rapa_kfwd, rapa_krev])
Rule('mTOR_Inhibitory_Complex', RAPA(s=2)%FKBP(sRAPA=2, sTOR=None) + mTOR(sRHEB=None, sATP=None, kFRB=None) | RAPA(s=2)%FKBP(sRAPA=2, sTOR=4)%mTOR(sRHEB=None, sATP=None, kFRB=4), *[rapa_comp_kfwd, rapa_comp_krev])
bind(mTOR, 'kFRB', RAPA, 's', [rapa_mt_kfwd, rapa_mt_krev])

#mTOR Binding ATP & Activated RHEB --> Activation Conditions w/ Binding of Both, RHEB 1st

Rule('mTOR_pRHEB_Binding', mTOR(sRHEB=None, sATP=None, kFRB=None) + pRHEB(p=None)| pRHEB(p=2)%mTOR(sRHEB=2, sATP=None, kFRB=None), *[kf_bind, kr_bind])
Rule('mTOR_pRHEB_to_ATP_Bind', pRHEB(p=2)%mTOR(sRHEB=2, sATP=None, kFRB=None) + ATP(s=None) | pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=None), *[kf_bind, kr_bind])

#Activated mTOR Binding Downstream S6K and 4EBP... CANT phosphorylate both S6K and 4EBP at the same time. Both occupy FRB domain

Rule('Activated_mTOR_S6K_Binding', pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=None) + S6K(s=None) | pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=4)%S6K(s=4), *[kf_bind, kr_bind])
Rule('Activated_mTOR_4EBP_Binding', pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=None) + EBP4(s=None) | pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=4)%EBP4(s=4), *[kf_bind, kr_bind])

#Production of mTOR Phosphorylated Prodcuts --> Results in Increased Transcription & Protein Translation 
Rule('Production_of_pS6K', pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=4)%S6K(s=4) >> RHEB(kAKT=None) + ADP(s=None) + mTOR(sRHEB=None, sATP=None, kFRB=None) + pS6K(p=None), *[kcat_s6])
Rule('Production_of_pEBP4', pRHEB(p=2)%ATP(s=3)%mTOR(sRHEB=2, sATP=3, kFRB=4)%EBP4(s=4) >> RHEB(kAKT=None) + ADP(s=None) + mTOR(sRHEB=None, sATP=None, kFRB=None) + pEBP4(p=None) , *[kcat_ebp4])

#pS6K Negative Feedback Loop --> Inhibits PI3K
Rule('pS6K_Binding_Unbound_PI3K', PI3K(sS6K=None, sIRS=None) + pS6K(p=None) | PI3K(sS6K=1, sIRS=None)%pS6K(p=1), *[kf_bind, kr_bind])
Rule('pS6K_Inactivating_Unbound_PI3K', PI3K(sS6K=1, sIRS=None)%pS6K(p=1) >> iPI3K(i=None) + S6K(s=None), *[kcat])


# Observables
Observable('opS6K', pS6K(p=None))
Observable('opEBP4', pEBP4(p=None))
Observable('opAKT', pAKT(p=None))
Observable('oINS', INS(s=None))
Observable('oFKBP', FKBP(sRAPA=None, sTOR=None))
Observable('oRAPA', RAPA(s=None))
Observable('oiPI3K', PI3K(sS6K=1, sIRS=None)%pS6K(p=1))
Observable('opPI3K', iPI3K(i=None))

#Simulate Pathway From 0 to 20 mins
t = np.linspace(0,1200)
out = odesolve(model, t)
#Print Observed Proteins
for i in range(0,len(out),5):
    print('pS6K: %.0f  pEBP4: %.0f  pAKT: %.0f  pPI3K: %.0f  FKBP: %.0f iPI3K: %.0f' % \
            (out['opS6K'][i],out['opEBP4'][i],out['opAKT'][i],out['opPI3K'][i],out['oFKBP'][i],out['oiPI3K'][i]))

# simres = ScipyOdeSimulator(model, tspan=t).run()
# yout = simres.all
# pl.ion()
# pl.figure()
# pl.plot(t, yout['opS6K'], label="pS6K")
# pl.plot(t, yout['opEBP4'], label="pEBP4")
# pl.plot(t, yout['oiPI3K'], label="inh_PI3K")
# pl.legend()
# pl.xlabel("Time (s)")
# pl.ylabel("Molecules/cell")
# pl.show()

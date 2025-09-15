# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 08:37:07 2024

@author: IGB
"""

import thermosteam as tmo
import biosteam as bst
from thermosteam import Chemical, functional as fn
from thermosteam.utils import chemical_cache

from biorefineries.HP.chemicals_data import HP_chemicals

#%% 
# Constants

_cal2joule=4.184

#%% 
@chemical_cache
def create_MeOH_chemicals():
    chems=tmo.Chemicals([])

    H2 = tmo.Chemical('H2')
    
    HCOOH = tmo.Chemical('HCOOH')

    C6H15N = tmo.Chemical('triethylamine')

    # C6H15N+HCOOH adduct modeling
    TREAHCOOH = tmo.Chemical.blank('TREAHCOOH', phase='l')
    TREAHCOOH.copy_models_from(C6H15N)
    TREAHCOOH.MW = 165.26
    TREAHCOOH.Tb = 800  # Set a fake high boiling point so that it never vaporizes during normal distillation or flash
    TREAHCOOH.Hf = -629e3  # Estimated
    TREAHCOOH.Psat.methods = []
    TREAHCOOH.Psat.tabulated = False

    # Prevent vaporization
    TREAHCOOH.Psat.methods = []
    TREAHCOOH.Psat.tabulated = False
    TREAHCOOH.Hvap.methods = []
    TREAHCOOH.Hvap.tabulated = False

    nBIM = tmo.Chemical('nBIM', search_ID='1-butylimidazole')

    # nBIM+HCOOH adduct modeling
    nBIMHCOOH = tmo.Chemical.blank('nBIMHCOOH', phase='l')
    nBIMHCOOH.copy_models_from(nBIM)
    nBIMHCOOH.MW = 170.22
    nBIMHCOOH.Tb = 800 # decomposes before boiling
    nBIMHCOOH.Pc = 3e6
    nBIMHCOOH.Hf = -625e3
    nBIMHCOOH.Psat.add_method(lambda T: 1e-10, name='const', Tmin=300, Tmax=800)

    phase_change_chemicals = ['HCOOH']

    for chem in chems:
        if chem.ID in phase_change_chemicals: pass
        elif chem.locked_state: pass
        else: 
            # Set phase_ref to avoid missing model errors
            if chem.phase_ref == 'g':
                chem.at_state('g')
            if chem.phase_ref == 'l':
                chem.at_state('l')
            if chem.phase_ref == 's':
                chem.at_state('s')
                
        
    
    chems.compile()
    return chems

#%% Add 3HP chemicals




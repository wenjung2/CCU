# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 15:12:18 2026

@author: IGB
"""


import CCU
import numpy as np
import biosteam as bst

#%% For sys_MeOH_water_electrolyzer_renewable system
system = CCU.create_full_system(ID='sys_MeOH_water_electrolyzer_renewable', water_electrolyzer=True)
system.flowsheet.BT.satisfy_system_electricity_demand=False

CCU.load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798, # Average 2023 https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=CCU.CFs['GWP_100']['Electricity'],
                                      electricity_price=CCU.price['Electricity'])
u = system.flowsheet.unit
s = system.flowsheet.stream

feedstock_ID = 'cornstover'
feedstock = s.cornstover
product_stream = s.ethanol

U101 = u.U101
@U101.add_specification(run=True, args=[83333.33333,0.2,0.3611,0.206,0.2005,0.041])
def set_feedstock_composition(feedstock_dry_flow, water_content, glucan_dry, xylan_dry, lignin_dry, ash_dry):
    cs = U101.ins[0]
    
    specified_components = ('Water', 'Glucan', 'Xylan', 'Lignin', 'Ash')
    other_components = [i.ID for i in cs.available_chemicals if i.ID not in specified_components]
    
    a = water_content
    b = glucan_dry
    c = xylan_dry
    d = lignin_dry
    e = ash_dry
    
    # new_water_mass = a
    dry_mass = 1-a
    new_glucan = b * (1-a)
    new_xylan = c * (1-a)
    new_lignin = d * (1-a)
    new_ash = e * (1-a)
    
    feedstock_flow = feedstock_dry_flow / (1-a)
    
    old_other_imass = {i:cs.imass[i] / sum([cs.imass[j] for j in other_components]) 
                       for i in other_components}
    new_other_mass_total = feedstock_flow * (1 - a - new_glucan - new_xylan - new_lignin - new_ash)
    new_other_mass = {i: new_other_mass_total * old_other_imass[i]
                       for i in other_components}
    
    # scaled by new flowrate
    updated_flows = {'Water': feedstock_flow * a,
                     'Glucan': feedstock_flow * new_glucan,
                     'Xylan': feedstock_flow * new_xylan,
                     'Lignin': feedstock_flow * new_lignin,
                     'Ash': feedstock_flow * new_ash,
                     **new_other_mass
                     }
                            
    cs.reset_flow(units='kg/hr', **updated_flows)

# =============================================================================
# create TEA
# =============================================================================
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185

tea = CCU.CellulosicIncentivesTEA(system=system, IRR=0.10, duration=(2023, 2053),
               depreciation='MACRS7', income_tax=0.35, operating_days=0.9*365,
               lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
               startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
               startup_VOCfrac=0.75, WC_over_FCI=0.05, finance_interest=0.08,
               finance_years=10, finance_fraction=0.6, 
               OSBL_units=(u.BT, u.CT, u.CWP, u.PWC, u.ADP, u.FWT, u.CIP),
               warehouse=0.04, site_development=0.09, additional_piping=0.045,
               proratable_costs=0.10, field_expenses=0.10, construction=0.20,
               contingency=0.10, other_indirect_costs=0.10, 
               labor_cost=3651112*get_flow_tpd()/2205,
               labor_burden=0.9, property_insurance=0.007, maintenance=0.03,
               steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT,
               carbon_credit=85, credit_years=12,)
               
system.operating_hours = tea.operating_days * 24
get_annual_factor = lambda: tea.operating_days * 24

feedstock.price = 0.091146891
s.sulfuric_acid.price = 0.145967079
s.ammonia.price = 0.875802473
s.cellulase.price = 0.408742354
s.DAP.price = 1.200932679
s.CSL.price = 0.099133292
s.caustic.price = 0.912957729
s.denaturant.price = 0.919936623
s.cooling_tower_chemicals.price = 2.144704919
s.FGD_lime.price = 0.13514
s.boiler_chemicals.price = 3.563518793
s.makeup_process_water.price = 0.001092283
s.makeup_RO_water.price = 0.002265475
u.BT.ash_disposal_price = -0.055494448
u.M301.solids_loading = 0.2
u.M301.enzyme_loading = 0.02
u.R303.saccharification_split = 0.1
u.R303.saccharification[2].X = 0.9
u.R303.cofermentation[0].X = 0.95
u.R303.cofermentation[4].X = 0.85
u.R302.glucose_to_ethanol.X = 0.9
u.R302.xylose_to_ethanol.X = 0.8
u.BT.boiler_efficiency = 0.8
u.BT.turbogenerator_efficiency = 0.85
s.natural_gas.price = 0.289367356
s.makeup_MEA.price = 1.44
s.catalyst_MeOH.price = 32.48578024


s.O2.price = 0.23

# s.MeOH.price = 0.8


# Contour plotting
y_data = np.linspace(85, 200, 9)
x_data = np.linspace(0.3,0.8,6)
w_data = []

def MSP_at_x_and_y_1(x,y):
    tea.carbon_credit = y
    s.MeOH.price = x
    for i in range(3): system.simulate()
    MSP = tea.solve_price(product_stream) / (product_stream.imass['Ethanol'] / product_stream.F_mass)
    return MSP

results = []
for j in y_data:
    w_data.append([])
    for i in x_data:
        try:
            print(i, j, MSP_at_x_and_y_1(i,j))
            w_data[-1].append(MSP_at_x_and_y_1(i,j))        
        except:
            print('Needs_interpolation')
            w_data[-1].append(0)  




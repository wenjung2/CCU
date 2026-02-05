# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 15:12:18 2026

@author: IGB
"""


import CCU

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
@U101.add_specification(run=True, args=[83333.33333,0.2,0.349,0.227,0.207,0.045])
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


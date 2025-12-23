# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 11:09:04 2025

@author: IGB
"""


import os
import numpy as np
import pandas as pd
import chaospy as cp
import biosteam as bst
from biosteam.process_tools import UnitGroup
from biosteam.evaluation import Model, Metric
from biosteam.evaluation._model import create_function
import CCU
from warnings import filterwarnings; filterwarnings('ignore')

available_systems = ['sys_ethanol_conventional', 'sys_ethanol_renewable', # in terms of electricity
                     'sys_MeOH_water_electrolyzer_renewable', 'sys_MeOH_water_electrolyzer_NG',
                     'sys_MeOH_hydrogen_conventional', 'sys_MeOH_hydrogen_renewable']

system_element_mapping = {available_systems[0]: {'A', 'D1'},
                          available_systems[1]: {'A', 'D2'},
                          available_systems[2]: {'A', 'B', 'C1', 'D2'},
                          available_systems[3]: {'A', 'B', 'C1', 'D1'},
                          available_systems[4]: {'A', 'B', 'C2', 'D1'},
                          available_systems[5]: {'A', 'B', 'C3', 'D1'},
                          }
fixed_parameter = 'Electricity unit price (renewable)'
fixed_value = 0.02

#%%
def create_model(system_name):
    if system_name == available_systems[0]:
        system = CCU.create_ethanol_system(ID='sys_ethanol_conventional')
    elif system_name == available_systems[1]:
        system = CCU.create_ethanol_system(ID='sys_ethanol_renewable')
    elif system_name == available_systems[2]:
        system = CCU.create_full_system(ID='sys_MeOH_water_electrolyzer_renewable', water_electrolyzer=True)
        system.flowsheet.BT.satisfy_system_electricity_demand=False
        # @system.flowsheet.PWC.add_specification(run=True)
        # def update_water_streams():
        #     u = system.flowsheet.unit
        #     s = system.flowsheet.stream
        #     u.PWC.makeup_water_streams = (u.CT.ins[1], u.BT.ins[2])
        #     u.PWC.process_water_streams = (s.warm_process_water_1, s.ammonia_process_water,\
        #                                    s.pretreatment_steam, s.warm_process_water_2,\
        #                                        s.saccharification_water, s.stripping_water,\
        #                                            u.S401.ins[1], u.R1101.ins[0], u.U1301.ins[2],\
        #                                                u.CIP.ins[0], u.FWT.ins[0])
    elif system_name == available_systems[3]:
        system = CCU.create_full_system(ID='sys_MeOH_water_electrolyzer_NG', water_electrolyzer=True)
        system.flowsheet.BT.satisfy_system_electricity_demand=True
        # @system.flowsheet.PWC.add_specification(run=True)
        # def update_water_streams():
        #     u = system.flowsheet.unit
        #     s = system.flowsheet.stream
        #     u.PWC.makeup_water_streams = (u.CT.ins[1], u.BT.ins[2])
        #     u.PWC.process_water_streams = (s.warm_process_water_1, s.ammonia_process_water,\
        #                                    s.pretreatment_steam, s.warm_process_water_2,\
        #                                        s.saccharification_water, s.stripping_water,\
        #                                            u.S401.ins[1], u.R1101.ins[0], u.U1301.ins[2],\
        #                                                u.CIP.ins[0], u.FWT.ins[0])
    elif system_name == available_systems[4]:
        system = CCU.system_hydrogen_purchased(ID='sys_MeOH_hydrogen_conventional', water_electrolyzer=False)
    elif system_name == available_systems[5]:
        system = CCU.system_hydrogen_purchased(ID='sys_MeOH_hydrogen_renewable', water_electrolyzer=False)
    
    
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

    # Unit groups (ethanol production + CCU + OSBL)
    ethanol_production = UnitGroup('ethanol_production', units = [i for i in system.units if \
                                                                  len(i.ID) < 5 and \
                                                                      i.ID[1] in ('1', '2', '3', '4') \
                                                                          and i.ID not in ('M1', 'M2')])
                                                             
    methanol_production = UnitGroup('methanol_production', units = [i for i in system.units if \
                                                                    len(i.ID) > 4 and (                 
                                                                    (i.ID[1] == '1' and i.ID[2] == '1') or \
                                                                    (i.ID[1] == '1' and i.ID[2] == '3'))])
                                              
    WWT = UnitGroup('WWT', units = [i for i in system.units if i.ID[1] == '6' \
                                        or i.ID in ('M1', 'M2', 'WWTC')])

    BT = UnitGroup('BT', units = (u.BT,)) 

    CT = UnitGroup('CT', units = (u.CT,))

    CWP = UnitGroup('CWP', units = (u.CWP,))

    PWC = UnitGroup('PWC', units = (u.PWC,))

    other_OSBL = UnitGroup('other_OSBL', units = [i for i in system.units if \
                                                      i not in ethanol_production \
                                                          and i not in methanol_production \
                                                              and i not in WWT \
                                                                  and i.ID not in ('BT', 'CT', 'CWP', 'PWC')])

    process_groups = [ethanol_production, methanol_production, \
                      WWT, BT, CT, CWP, PWC, other_OSBL]

    process_groups_dict = {}
    for i in range(len(process_groups)):
        group = process_groups[i]
        process_groups_dict[group.name] = group
    
    # =============================================================================
    # add system specification capable of setting feedstock components
    # =============================================================================
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

    # tea = CCU.EtOH_TEA(system=system, IRR=0.10, duration=(2023, 2053),
    #                depreciation='MACRS7', income_tax=0.35, operating_days=0.9*365,
    #                lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
    #                startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
    #                startup_VOCfrac=0.75, WC_over_FCI=0.05, finance_interest=0.08,
    #                finance_years=10, finance_fraction=0.6, 
    #                OSBL_units=(u.BT, u.CT, u.CWP, u.PWC, u.ADP, u.FWT, u.CIP),
    #                warehouse=0.04, site_development=0.09, additional_piping=0.045,
    #                proratable_costs=0.10, field_expenses=0.10, construction=0.20,
    #                contingency=0.10, other_indirect_costs=0.10, 
    #                labor_cost=3651112*get_flow_tpd()/2205,
    #                labor_burden=0.9, property_insurance=0.007, maintenance=0.03,
    #                steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)
    
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

    # =============================================================================
    # create LCA
    # =============================================================================
                                  
    if system_name == available_systems[0] or system_name == available_systems[1]:
        CCU.CFs['GWP_100']['O2'] = 0 # set once from dic, not function not param
        input_biogenic_carbon_streams = (feedstock, s.cellulase, s.CSL,
                                         s.natural_gas)
        by_products = []  # No coproducts for ethanol system
    elif system_name == available_systems[2] or system_name == available_systems[3]:
        input_biogenic_carbon_streams = (feedstock, s.cellulase, s.CSL, s.makeup_MEA,
                                         s.natural_gas)
        by_products = [s.MeOH, s.O2]  # MeOH and O2
    else:
        CCU.CFs['GWP_100']['O2'] = 0
        input_biogenic_carbon_streams = (feedstock, s.cellulase, s.CSL, s.makeup_MEA,
                                         s.natural_gas)
        by_products = [s.MeOH] 

    lca = CCU.create_CCU_lca(system=system,
                             CFs=CCU.CFs,
                             feedstock=feedstock,
                             feedstock_ID=feedstock_ID,
                             main_product=product_stream,
                             main_product_chemical_IDs=['Ethanol'],
                             by_products=by_products,
                             cooling_tower=u.CT,
                             chilled_water_processing_units=u.CWP,
                             boiler=u.BT, has_turbogenerator=True,
                             add_EOL_GWP=True,
                             input_biogenic_carbon_streams=input_biogenic_carbon_streams)
                         

    # =============================================================================
    # Metrics
    # ============================================================================
    # 0. TEA
    def get_MSP():
        for i in range(3):
            product_price = tea.solve_price(product_stream)
        return product_price

    get_yield = lambda: product_stream.F_mass * get_annual_factor() / 1e6 # in 1000 MT
    get_purity = lambda: product_stream.imass['Ethanol'] / product_stream.F_mass

    # adjusted for purity
    get_adjusted_MSP = lambda: get_MSP() / get_purity()
    get_adjusted_yield = lambda: get_yield() / get_purity()

    get_overall_TCI = lambda: tea.TCI / 1e6
    get_overall_installed_cost = lambda: tea.installed_equipment_cost / 1e6

    # annual operating cost
    get_overall_AOC = lambda: tea.AOC / 1e6
    get_material_cost = lambda: (tea.material_cost + abs(s.ash.F_mass * s.ash.price)) / 1e6
    get_overall_FOC = lambda: tea.FOC / 1e6
    
    get_hydrogen_AOC = lambda: s.hydrogen.cost * get_annual_factor() / 1e6
    get_NG_AOC = lambda: s.natural_gas.cost * get_annual_factor() / 1e6
    # annual sale revenue from products, note that electricity credit is not included,
    # but negative sales from waste disposal are included
    # (i.e., wastes are products of negative selling price)
    get_annual_sale = lambda: tea.sales / 1e6

    # system power usage 
    get_excess_electricity = lambda: system.get_electricity_production() - system.get_electricity_consumption() # kWh per year
    get_electricity_revenue = lambda: get_excess_electricity() * bst.PowerUtility.price / 1e6 # 10^6 $ per year
    

    metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'TEA'),
               Metric('Production rate', get_yield, '10^6 kg/yr', 'TEA'),
               Metric('Product purity', get_purity, '%', 'TEA'),
               Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg', 'TEA'),
               Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr', 'TEA'),
               Metric('Total capital investment', get_overall_TCI, '10^6 $', 'TEA'),
               Metric('Total installed equipment cost', get_overall_installed_cost, '10^6 $', 'TEA'),
               Metric('Annual operating cost (incl. electricity credit)', get_overall_AOC, '10^6 $/yr', 'TEA'),
               Metric('Annual material cost (incl. boiler ash disposal)', get_material_cost, '10^6 $/yr', 'TEA'),
               Metric('Fixed operating cost', get_overall_FOC, '10^6 $/yr', 'TEA'),
               Metric('NG cost', get_NG_AOC, '10^6 $/yr', 'TEA'),
               Metric('Annual product sale (excl. electricity)', get_annual_sale, '10^6 $/yr', 'TEA'),
               Metric('Annual electricity output', get_excess_electricity, 'kWh/yr', 'TEA'),
               Metric('Annual electricity revenue', get_electricity_revenue, '10^6 $/yr', 'TEA'),]
    
               
                      
    
    for i in process_groups: i.autofill_metrics(shorthand=False, 
                                                electricity_production=False, 
                                                electricity_consumption=True,
                                                material_cost=True)
    metrics_labels_dict = {
        'Installed cost':(0, '10^6 $'), 
        'Material cost':(4,'USD/h'), 
        'Cooling duty':(1,'GJ/h'), 
        'Heating duty':(2,'GJ/h'), 
        'Electricity usage':(3, 'MW'),}
    
    for m, u_i in metrics_labels_dict.items():
        for ug in process_groups:
            metrics.append(Metric(ug.name, ug.metrics[u_i[0]], u_i[1], m))

    # 1. Carbon use
    all_products = [product_stream] + by_products
    emissions = [i for i in s if i.source and not i.sink and i not in all_products]
    
    check_C_balance = lambda: lca.system_carbon_balance
    
    get_C_emissions_per_ethanol = lambda: get_C_emissions() / get_C_ethanol()
    get_carbon_efficiency = lambda: get_C_ethanol() / get_C_feedstock()
    get_carbon_efficiency_total = lambda: get_C_ethanol() / get_C_in()
    
    get_C_in = lambda: sum([i.get_atomic_flow('C') for i in system.feeds])
    get_C_feedstock = lambda: feedstock.get_atomic_flow('C')
    get_C_CSL = lambda: system.flowsheet.CSL.get_atomic_flow('C')
    get_C_cellulase = lambda: system.flowsheet.cellulase.get_atomic_flow('C')
    get_C_denaturant = lambda: system.flowsheet.denaturant.get_atomic_flow('C')
    
    get_C_ethanol = lambda: system.flowsheet.ethanol.get_atomic_flow('C')
   
    get_C_emissions = lambda: sum([i.get_atomic_flow('C') for i in emissions])
    get_C_emissions_fermentation = lambda: system.flowsheet.D401.outs[0].get_atomic_flow('C')
    get_C_emissions_boiler = lambda: system.flowsheet.BT.outs[0].get_atomic_flow('C')
    get_C_emissions_WWT = lambda: system.flowsheet.R602.outs[0].get_atomic_flow('C') + system.flowsheet.S604.outs[1].get_atomic_flow('C')
    
    
    
    metrics.extend((Metric('Carbon balance', check_C_balance, '', 'Carbon'),))
    metrics.extend((Metric('C emissions per ethanol', get_C_emissions_per_ethanol, '', 'Carbon'),))
    metrics.extend((Metric('C efficiency', get_carbon_efficiency, '', 'Carbon'),))
    metrics.extend((Metric('Total C efficiency', get_carbon_efficiency_total, '', 'Carbon'),))
    metrics.extend((Metric('Total C in', get_C_in, '', 'Carbon'),))
    metrics.extend((Metric('Feedstock C', get_C_feedstock, '', 'Carbon'),))
    metrics.extend((Metric('CSL C', get_C_CSL, '', 'Carbon'),))
    metrics.extend((Metric('Cellulase C', get_C_cellulase, '', 'Carbon'),))
    metrics.extend((Metric('Denaturant C', get_C_denaturant, '', 'Carbon'),))
    metrics.extend((Metric('Ethanol C', get_C_ethanol, '', 'Carbon'),))
    metrics.extend((Metric('Total emissions C', get_C_emissions, '', 'Carbon'),))
    metrics.extend((Metric('Fermentation emissions C', get_C_emissions_fermentation, '', 'Carbon'),))
    metrics.extend((Metric('Boiler emissions C', get_C_emissions_boiler, '', 'Carbon'),))
    metrics.extend((Metric('WWT emissions C', get_C_emissions_WWT, '', 'Carbon'),))
    
    # metrics for CCU
    if system_name in [available_systems[2], available_systems[3], available_systems[4], available_systems[5]]:
        
        get_S1300_split = lambda: system.flowsheet.S1300.split[0]
        
        get_C_fermentation_used = lambda: system.flowsheet.M1302.ins[0].get_atomic_flow('C')
        get_C_boiler_used = lambda: system.flowsheet.M1302.ins[1].get_atomic_flow('C')
        get_C_natural_gas = lambda: system.flowsheet.BT.ins[3].get_atomic_flow('C')
        get_C_methanol = lambda: system.flowsheet.MeOH.get_atomic_flow('C')
        get_C_emissions_MeOH_offgas = lambda: system.flowsheet.gas_out.get_atomic_flow('C')
        
        get_C_emissions_per_ethanol_MeOH = lambda: get_C_emissions() / (get_C_ethanol() + get_C_methanol())
        
        get_carbon_efficiency_CCU = lambda:  (get_C_ethanol() + get_C_methanol()) / get_C_feedstock()
        get_carbon_efficiency_total_CCU = lambda:  (get_C_ethanol() + get_C_methanol()) / get_C_in()
        
        get_45Q_incentives = lambda: tea.annual_credit
        
        metrics.extend((Metric('S1300 split input', get_S1300_split, '', 'Carbon'),))
        metrics.extend((Metric('Fermentation C used', get_C_fermentation_used, '', 'Carbon'),))
        metrics.extend((Metric('Boiler C used', get_C_boiler_used, '', 'Carbon'),))
        metrics.extend((Metric('Natural gas C', get_C_natural_gas, '', 'Carbon'),))
        metrics.extend((Metric('MeOH C', get_C_methanol, '', 'Carbon'),))
        metrics.extend((Metric('MeOH_offgas C', get_C_emissions_MeOH_offgas, '', 'Carbon'),))
        metrics.extend((Metric('C emissions per ethanol MeOH', get_C_emissions_per_ethanol_MeOH, '', 'Carbon'),))
        metrics.extend((Metric('C efficiency CCU', get_carbon_efficiency_CCU, '', 'Carbon'),))
        metrics.extend((Metric('Total C efficiency CCU', get_carbon_efficiency_total_CCU, '', 'Carbon'),))
        metrics.extend((Metric('CCU incentives', get_45Q_incentives, '', 'Carbon'),))
    else:
        pass
    
    if system_name in [available_systems[2]]:
        # Changes in TEA
        get_electricity_input = lambda: (system.get_electricity_consumption() - system.get_electricity_production()) / get_annual_factor() # kWh per hour
        get_annual_electricity_cost = lambda: (system.get_electricity_consumption() - system.get_electricity_production()) * bst.PowerUtility.price / 1e6
        get_electricity_consumption_R1101 = lambda: system.flowsheet.R1101.power_utility.rate  # kWh per hour
        get_hydrogen_flow_R1101 = lambda: system.flowsheet.R1101.outs[0].imass['H2']
        get_normalized_hydrogen_power_consump = lambda: get_electricity_consumption_R1101() / get_hydrogen_flow_R1101()
        check_electricity = lambda: get_electricity_input() / get_electricity_consumption_R1101()
        
        metrics.extend((Metric('Electricity input', get_electricity_input, 'kWh/hr', 'TEA'),))
        metrics.extend((Metric('Electricity annual cost', get_annual_electricity_cost, 'MM$/year', 'TEA'),))
        metrics.extend((Metric('Electricity input for R1101', get_electricity_consumption_R1101, 'kWh/hr', 'TEA'),))
        metrics.extend((Metric('Electricity input per H2', get_normalized_hydrogen_power_consump, 'kWh/kg', 'TEA'),))
        metrics.extend((Metric('Electricity check', check_electricity, '', 'TEA'),))
    else:
        pass

    # 2. LCA
    if system_name == available_systems[0] or system_name == available_systems[1]:
        get_GWP = lambda: lca.GWP
        get_material_GWP = lambda: lca.material_GWP
        get_other_materials_GWP = lambda: get_material_GWP() - lca.material_GWP_breakdown['CSL'] -\
            lca.material_GWP_breakdown['DAP'] - lca.material_GWP_breakdown['CH4'] -\
                lca.material_GWP_breakdown['Cellulase']
        
        # using energy allocation
        sec_per_hr = 60 * 60
        kJ_per_GGE = 120276
        GGE_electricity = lambda: sec_per_hr/kJ_per_GGE * get_excess_electricity()
        GGE_ethanol = lambda: s.ethanol.get_property('LHV', 'GGE/hr') * system.operating_hours
        ethanol_energy_allocation = lambda: GGE_ethanol() / (GGE_electricity() + GGE_ethanol())
        electricity_energy_allocation = lambda: GGE_electricity() / (GGE_electricity() + GGE_ethanol())
        
        ethanol_GWP_by_energy = lambda: get_GWP_before_electricity_offset() * ethanol_energy_allocation()
        electricity_GWP_by_energy = lambda: get_GWP_before_electricity_offset() * electricity_energy_allocation()
        
        metrics.append(Metric('GWP100a - ethanol by allocation', ethanol_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - electricity by allocation', electricity_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
        
    elif system_name == available_systems[2] or system_name == available_systems[3]:
        get_GWP = lambda: lca.GWP - lca.material_GWP_breakdown['O2']
        get_material_GWP = get_material_GWP_no_O2 = lambda: lca.material_GWP - lca.material_GWP_breakdown['O2']
        get_other_materials_GWP = lambda: get_material_GWP_no_O2() -\
            lca.material_GWP_breakdown['CSL'] - lca.material_GWP_breakdown['DAP'] - lca.material_GWP_breakdown['CH4'] -\
                lca.material_GWP_breakdown['MEA'] - lca.material_GWP_breakdown['Cellulase']
        Metric('H2 cost', get_hydrogen_AOC, '10^6 $/yr', 'TEA'),
        metrics.append(Metric('GWP100a - Coproduct credit - Methanol', lambda: lca.GWP_byproduct_credit(0), 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - Coproduct credit - O2', lambda: lca.GWP_byproduct_credit(1), 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - Coproduct credit - total', lambda: lca.GWP_byproduct_credit_total(), 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP - O2', lambda: CCU.CFs['GWP_100']['O2'], 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('Amount - O2', lambda: s.O2.imass['O2'], 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('Amount - ETOH', lambda: s.ethanol.imass['Ethanol'], 'kg-CO2-eq/kg', 'LCA'))
        
        # using hybrid allocation (O2 displaced, MeOH and EtOH energy allocation)
        GWP_without_EtOH = lambda: get_GWP() + lca.GWP_byproduct_credit(0)
        ethanol_energy_allocation = lambda: s.ethanol.get_property('LHV', 'GGE/hr') / (s.ethanol.get_property('LHV', 'GGE/hr')
                                                                                       +s.MeOH.get_property('LHV', 'GGE/hr'))
        MeOH_energy_allocation = lambda: s.MeOH.get_property('LHV', 'GGE/hr') / (s.ethanol.get_property('LHV', 'GGE/hr')
                                                                                       +s.MeOH.get_property('LHV', 'GGE/hr'))
        ethanol_GWP_by_energy = lambda: GWP_without_EtOH() * ethanol_energy_allocation()
        MeOH_GWP_by_energy = lambda: GWP_without_EtOH() * MeOH_energy_allocation()
        
        metrics.append(Metric('GWP100a - ethanol by allocation', ethanol_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - MeOH by allocation', MeOH_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
    else:
        get_GWP = lambda: lca.GWP
        get_material_GWP = lambda: lca.material_GWP
        get_other_materials_GWP = lambda: get_material_GWP() - lca.material_GWP_breakdown['CSL'] -\
            lca.material_GWP_breakdown['DAP'] - lca.material_GWP_breakdown['CH4'] -\
                lca.material_GWP_breakdown['Cellulase'] - lca.material_GWP_breakdown['H2']
        Metric('H2 cost', get_hydrogen_AOC, '10^6 $/yr', 'TEA'),
        metrics.append(Metric('GWP100a - Coproduct credit - Methanol', lambda: lca.GWP_byproduct_credit(0), 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - Materials breakdown - H2', lambda: lca.material_GWP_breakdown['H2'], 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP - H2', lambda: CCU.CFs['GWP_100']['H2'], 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('Amount - H2', lambda: s.hydrogen.imass['H2'], 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('Amount - ETOH', lambda: s.ethanol.imass['Ethanol'], 'kg-CO2-eq/kg', 'LCA'))
        
        # using energy allocation (MeOH and EtOH)
        ethanol_energy_allocation = lambda: s.ethanol.get_property('LHV', 'GGE/hr') / (s.ethanol.get_property('LHV', 'GGE/hr')
                                                                                       +s.MeOH.get_property('LHV', 'GGE/hr'))
        MeOH_energy_allocation = lambda: s.MeOH.get_property('LHV', 'GGE/hr') / (s.ethanol.get_property('LHV', 'GGE/hr')
                                                                                       +s.MeOH.get_property('LHV', 'GGE/hr'))
        ethanol_GWP_by_energy = lambda: get_GWP() * ethanol_energy_allocation()
        MeOH_GWP_by_energy = lambda: get_GWP() * MeOH_energy_allocation()
        
        metrics.append(Metric('GWP100a - ethanol by allocation', ethanol_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
        metrics.append(Metric('GWP100a - MeOH by allocation', MeOH_GWP_by_energy, 'kg-CO2-eq/kg', 'LCA'))
        
    get_GWP_before_electricity_offset = lambda: get_GWP() - lca.net_electricity_GWP
    
    metrics.append(Metric('Total GWP100a', get_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('Total GWP100a before electricity offset', get_GWP_before_electricity_offset, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Electricity', lambda: lca.net_electricity_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Direct emissions', lambda: lca.direct_emissions_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Direct biogenic emissions', lambda: lca.biogenic_emissions_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - EoL emissions', lambda: lca.EOL_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Direct non-biogenic emissions', lambda: lca.direct_non_biogenic_emissions_GWP, 'kg-CO2-eq/kg', 'LCA'))

    metrics.append(Metric('GWP100a - Feedstock (FGHTP)', lambda: lca.FGHTP_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Materials (except feedstock)', get_material_GWP, 'kg-CO2-eq/kg', 'LCA'))
    
    metrics.append(Metric('GWP100a - Materials other', get_other_materials_GWP, 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Materials breakdown - CSL', lambda: lca.material_GWP_breakdown['CSL'], 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Materials breakdown - DAP', lambda: lca.material_GWP_breakdown['DAP'], 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Materials breakdown - CH4', lambda: lca.material_GWP_breakdown['CH4'], 'kg-CO2-eq/kg', 'LCA'))
    metrics.append(Metric('GWP100a - Materials breakdown - Cellulase', lambda: lca.material_GWP_breakdown['Cellulase'], 'kg-CO2-eq/kg', 'LCA'))
    
    # =============================================================================
    # Set up model
    # =============================================================================
    
    model = Model(system, metrics)
    param = model.parameter
    
    elements_for_system = system_element_mapping.get(system_name, set())
    
    # Generate the required namespace
    namespace_dict = {}
    exclude_from_globals = [
        'search',
        'register',
        'register_safely',
        'discard',
        'clear',
        'mark_safe_to_replace',
        'unmark_safe_to_replace']

    namespace_dict.update({k:s.__getitem__(k) for k in s.__dir__() if not k in exclude_from_globals})
    namespace_dict.update({k:u.__getitem__(k) for k in u.__dir__() if not k in exclude_from_globals})
    namespace_dict['feedstock'] = feedstock
    
    namespace_dict['tea'] = tea
    namespace_dict['lca'] = lca

    PowerUtility = bst.PowerUtility
    namespace_dict['PowerUtility'] = PowerUtility
    
    namespace = system.flowsheet.to_dict() | namespace_dict
    
    for i, row in dist_table.iterrows():
        param_name = row['Parameter name']
        element = row['Element']
        kind = row['Kind']
        units = row['Units']
        baseline = row['Baseline']
        statement = row['Statement']

        if element not in elements_for_system:
            continue
            

        param(name=param_name,
              setter=create_function(statement, namespace), 
              element=element, 
              kind=kind, 
              units=units,
              baseline=baseline,
              distribution=None)
    # =============================================================================
    # Bugfix barrage
    # =============================================================================
    def reset_and_reload():
        print('Resetting cache and emptying recycles ...')
        system.reset_cache()
        system.empty_recycles()
    def reset_and_switch_solver(solver_ID):
        system.reset_cache()
        system.empty_recycles()
        system.converge_method = solver_ID
        print(f"Trying {solver_ID} ...")
    def run_bugfix_barrage():
        try:
            reset_and_reload()
            system.simulate()
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('fixedpoint')
                system.simulate()
            except Exception as e:
                print(str(e))
                try:
                    reset_and_switch_solver('aitken')
                    system.simulate()
                except Exception as e:
                    print(str(e))
                    # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                    print("Bugfix barrage failed.\n")
                    # breakpoint()
                    raise e

    def model_specification():
        try:
            system.simulate()
        except Exception as e:
            str_e = str(e).lower()
            print('Error in model spec: %s'%str_e)
            run_bugfix_barrage()
    model.specification = model_specification()
    return model

#%% Generate parameters and samples

if __name__ == '__main__':
    EtOH_MeOH_filepath = CCU.EtOH.__file__.replace('\\__init__.py', '')
    input_folder = os.path.join(EtOH_MeOH_filepath, 'analyses', 'parameter_distributions')
    os.makedirs(input_folder, exist_ok=True)
    
    full_parameter_distributions_filename = 'parameter-distributions_full.xlsx'
    full_parameter_distributions_filepath = os.path.join(input_folder, full_parameter_distributions_filename)
    
    # Read the table
    dist_table  = pd.read_excel(full_parameter_distributions_filepath) 
    
    param_names = dist_table['Parameter name'].tolist()
    elements = dist_table['Element'].tolist()
    
    distributions = []
    sample_positions = []
    
    for pos, (_, row) in enumerate(dist_table.iterrows()):
        name = str(row["Parameter name"]).strip()
        
        if name == fixed_parameter:
            continue
        
        shape = row['Shape'].strip().lower()
        lower = row['Lower']
        upper = row['Upper']
    
        if shape == 'triangular':
            mode = row['Midpoint']
            dist = cp.Triangle(lower, mode, upper)
        elif shape == 'uniform':
            dist = cp.Uniform(lower, upper)
        else:
            raise ValueError(f"Unsupported shape: {shape}")
    
        distributions.append(dist)
        sample_positions.append(pos)
    
    # Generate N Latin Hypercube samples
    N = 10
    joint_dist = cp.J(*distributions)
    
    rand_samples = joint_dist.sample(size=N, rule="L", seed=3221).T 
    
    full_samples = np.empty((N, len(dist_table)), dtype=float)
    
    fixed_pos = param_names.index(fixed_parameter)
    full_samples[:, fixed_pos] = fixed_value
    for j, pos in enumerate(sample_positions):
        full_samples[:, pos] = rand_samples[:, j]
    sample_df = pd.DataFrame(full_samples)
    # Shift samples one column right by adding empty column at index 0
    sample_df.insert(0, 'Empty', '')
    
    # Create column A data for sample indices (rows 2 onward)
    sample_indices = list(range(1, N + 1))
    
    # Prepare metadata rows (row 0 and 1) with empty first cell for alignment
    row1 = ['Element'] + elements
    row2 = ['Parameter name'] + param_names
    
    # Build full data as list of lists
    output_data = [row1, row2]
    
    # Append samples with index in column A
    for i, sample_row in enumerate(sample_df.values, start=1):
        # replace first cell (empty) with index i
        sample_row[0] = i
        output_data.append(sample_row.tolist())
    
    # Convert to DataFrame and export
    output_df = pd.DataFrame(output_data)
    output_filename = f'{N}_full_samples.xlsx'
    output_path = os.path.join(input_folder, output_filename)
    output_df.to_excel(output_path, index=False, header=False)
    
    #%%
    from datetime import datetime
    
    EtOH_MeOH_results_filepath = EtOH_MeOH_filepath + '\\analyses\\results\\'
    
    dateTimeObj = datetime.now()
    
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    
    def run_model(sys_name, notify_runs=10):
        model = create_model(system_name=sys_name)
        
        allowed_elements = system_element_mapping.get(sys_name, set())
    
        # Filter parameter names based on element
        allowed_param_names = [
            name for name, elem in zip(param_names, elements)
            if elem in allowed_elements
        ]
    
        sample_columns = ['Parameter name'] + param_names
        param_col_indices = [
            sample_columns.index(pname)
            for pname in allowed_param_names
        ]
    
        # Sample rows (skip metadata rows 0 and 1)
        sample_list = []
        for i in range(2, N + 2):  # row index in output_df
            row = output_df.iloc[i]
            sample_row = [row[j] for j in param_col_indices]
            sample_list.append(sample_row)
            
        sample_array = np.array(sample_list)
        
        model.load_samples(sample_array)
        
        # Baseline results
        baseline_initial = model.metrics_at_baseline()
        baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                                columns=baseline_initial.keys())

        model.evaluate(notify=notify_runs)
        
        # Percentiles
        percentiles = [0.05, 0.25, 0.50, 0.75, 0.95]
        percentiles_df = model.table.quantile(q=percentiles)
        
        # Spearman's rank correlation
        df_rho, df_p = model.spearman_r(filter='omit nan')
        
        target_indicators = ["Minimum selling price",
                             "Total gwp100a",]
        
        sig_params = CCU.get_significant_params(rho_df=df_rho, p_df = df_p, 
                                                indicator_filter = target_indicators)
            

        file_to_save = EtOH_MeOH_results_filepath\
            +'_' + sys_name + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
            + '_' + '_' + str(N) + 'sims'
        # Output to Excel
        with pd.ExcelWriter(file_to_save+'_'+'_1_full_evaluation.xlsx') as writer:
            baseline.to_excel(writer, sheet_name='Baseline')
            percentiles_df.to_excel(writer, sheet_name='Percentile results')
            sig_params.to_excel(writer, sheet_name='Significant parameters')
            df_rho.to_excel(writer, sheet_name='df_rho')
            df_rho.to_excel(writer, sheet_name='df_p')
            model.table.to_excel(writer, sheet_name='Raw data')
        
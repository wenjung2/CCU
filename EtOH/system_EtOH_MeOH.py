# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:38:17 2025

@author: IGB
"""
# this runs in succinic environment.
    # 1. conda activate ccu
    # 2. Select file path (biosteam/thermosteam) here
    
import biosteam as bst
import CCU
from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system

#%% Captured CO2 to MeOH

@bst.SystemFactory(ID='sys_ccu',
                   ins=[dict(ID='emissions_BT', CO2=1),
                        dict(ID='emissions_fermentation', CO2=1),
                        dict(ID='makeup_MEA', H2O=1),
                        dict(ID='makeup_water_1', H2O=1),
                        dict(ID='hydrogen', H2=1),
                        dict(ID='water_stream_1', H2O=1),
                        dict(ID='catalyst_MeOH', CaSO4=1),],
                   outs=[dict(ID='bottom_water', H2O=1),
                         dict(ID='oxygen', O2=1),
                         dict(ID='spent_catalyst', CaSO4=1),
                         dict(ID='gas_out', CO2=1),
                         dict(ID='MeOH', CH3OH=1),])
def create_ccu_system(ins, outs, water_electrolyzer=None, hydrogen_green=None, hydrogen_blue=None, hydrogen_gray=None):
    emissions_BT, emissions_fermentation, makeup_MEA, makeup_water_1, hydrogen, water_stream_1, catalyst_MeOH = ins
    bottom_water, oxygen, spent_catalyst, gas_out, MeOH = outs
    
    # capture CO2
    S1300 = bst.Splitter('S1300', ins=emissions_BT, outs=('captured', 'purge'), split=1)
    S1300.prioritize = True
    S1300.register_alias('carbon_capture_splitter')
    S1300.maximum_biorefinery_percent_capture = 150 # between 1 to <2
    @S1300.add_specification(run=True)
    def maximum_biorefinery_percent_capture(): # Percent of cellulosic biorefinery carbon emissions
        emissions = S1300.ins[0]
        try:
            baseline_biorefinery_emissions = S1300.baseline_biorefinery_emissions
        except:
            S1300.baseline_biorefinery_emissions = baseline_biorefinery_emissions = emissions.F_mol
            S1300.split[:] = 1
            return
        new_emissions = emissions.F_mol
        maximum_emissions = S1300.maximum_biorefinery_percent_capture * baseline_biorefinery_emissions / 100
        if new_emissions <= maximum_emissions:
            S1300.split[:] = 1
        else:
            S1300.split[:] = maximum_emissions / new_emissions
        # print(f"{maximum_emissions / new_emissions:.0%}")
    
    U1301 = bst.AmineAbsorption('U1301', ins=(S1300-0, makeup_MEA, makeup_water_1),\
                                  outs=('vent', 'concentrated'), CO2_recovery=0.8)
    U1301.outs[1].phase = 'g'
    M1302 = bst.Mixer('M1302', ins=(emissions_fermentation, U1301-1), outs='')
    S1301 = bst.Splitter('S1301', ins=M1302-0, outs=('concentrated_CO2', 'other_gases'), split=dict(CO2=1.0))
    
    ks = [bst.IsentropicCompressor(P=3*101325), 
          bst.IsentropicCompressor(P=8.8*101325),
          bst.IsentropicCompressor(P=26.3*101325)]

    hxs = [bst.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]
    
    # First 3-stage compression of CO2
    # C1101 = bst.MultistageCompressor('C1101', ins=S1301-0, outs='', compressors=ks, hxs=hxs) # won't work (no ins connection)
    C1101 = bst.MultistageCompressor('C1101', ins=S1301-0, outs='', n_stages=3, pr=2.9)
    
    # Final compression of CO2
    C1102 = bst.IsentropicCompressor('C1102', ins=C1101-0, outs='', P=78*101325)
    
    if water_electrolyzer:
        ins.remove(hydrogen)
        R1101 = CCU.Electrolyzer('R1101', ins=water_stream_1, outs=('', oxygen))
        @R1101.add_specification(run=True)
        def adjust_water_flow():
            U1301.run()
            R1101.ins[0].imol['H2O'] = U1301.outs[1].imol['CO2'] * 3 # H2:CO2 = 3:1
        # H2 compressed to same pressure
        C1103 = bst.IsentropicCompressor('C1103', ins=R1101-0, outs='', P=78*101325, vle=True)
    else:
        ins.remove(water_stream_1)
        outs.remove(oxygen)
        C1103 = bst.IsentropicCompressor('C1103', ins=hydrogen, outs='', P=78*101325, vle=True)
        @C1103.add_specification(run=True)
        def adjust_hydrogen_flow():
            U1301.run()
            C1103.ins[0].imol['H2'] = U1301.outs[0].imol['CO2'] * 3 # H2:CO2 = 3:1
        
    M1101 = bst.Mixer('M1101', ins=(C1102-0, C1103-0), outs='', rigorous=True)

    M1102 = bst.Mixer('M1102', ins=(M1101-0, ''), outs='')
    M1102.prioritize=True

    H1101 = bst.HXprocess('H1101', ins=(M1102-0, ''), outs=('', ''),)

    R1102 = CCU.MeOH_SynthesisReactor('R1102', ins=(H1101-0, catalyst_MeOH),\
                                         outs=('product', spent_catalyst))

    S1101 = bst.Splitter('S1101', ins=R1102-0, outs=(1-H1101, ''), split=0.6)

    # Couldn't model heat exchange process between S1101-1 and reboiler stream in DT1REB of paper 
    # and use a simple HX for T required in H1103
    H1101_1 = bst.HXutility('H1101_1', ins=S1101-1, outs='', T=156+273.15)

    H1103 = bst.HXprocess('H1103', ins=(H1101_1-0, ''), outs=('', ''), phase0='g', phase1='l')

    D1101 = bst.BinaryDistillation('D1101', ins=H1103-1, outs=('gas_MEOH', bottom_water),
                                     LHK=('CH3OH', 'H2O'),
                                     Lr=0.9999, Hr=0.9999, k=2,
                                     is_divided=True)
                                 
    C1105 = bst.IsentropicCompressor('C1105', ins=D1101-0, outs='', P=1.2*101325)

    H1104 = bst.HXutility('H1104', ins=C1105-0, outs='', T=40+273.15, rigorous=True)

    S1104 = bst.PhaseSplitter('S1104', ins=H1104-0, outs=('gas_final', ''))

    V1101 = bst.IsenthalpicValve('V1101', ins=H1101-1, outs='', P=73.6*101325)

    M1103 = bst.Mixer('M1103', ins=(V1101-0, H1103-0), outs='', rigorous=True)

    H1102 = bst.HXutility('H1102', ins=M1103-0, outs='', T=35+273.15, rigorous=True)

    S1102 = bst.PhaseSplitter('S1102', ins=H1102-0, outs=('gas', 'condensed_water_and_methanol'))

    S1103 = bst.Splitter('S1103', ins=S1102-0, outs=('purge', 'recycled'), split=0.01)

    C1104 = bst.IsentropicCompressor('C1104', ins=S1103-1, outs=1-M1102, P=78*101325)

    # For condensed water and methanol in S1102 (composed of methanol, water and residual dissolved gases)
    V1102 = bst.IsenthalpicValve('V1102', ins=S1102-1, outs='', P=10*101325)

    V1103 = bst.IsenthalpicValve('V1103', ins=V1102-0, outs='', P=1.2*101325)

    F1101 = bst.SplitFlash('F1101', ins=V1103-0, outs=('gas_F1101', 1-H1103), T=22+273.15, P=1.2*101325, split=dict(CO2=1.0,
                                                                                                           H2=1.0))

    M1104 = bst.Mixer('M1104', ins=(S1103-0, F1101-0, S1104-0), outs=gas_out)
    
    T1101 = bst.StorageTank('T1101', ins=S1104-1, outs=MeOH)
    
    
def create_full_system():
    chems = CCU.create_MeOH_chemicals()
    bst.settings.set_thermo(chems, cache=True)
    biorefinery_sys = CCU.create_cellulosic_ethanol_system('sys_ethanol_cs')
    
    emissions_BT = biorefinery_sys.flowsheet.BT.outs[0] # Boiler emissions
    emissions_fermentation = biorefinery_sys.flowsheet.D401.outs[0] # Fermentation emissions
    carbon_capture_sys = CCU.create_ccu_system(ins=[emissions_BT,
                                                emissions_fermentation,
                                                'makeup_MEA',
                                                'makeup_water_1',
                                                'hydrogen',
                                                'water_stream_1',
                                                'catalyst_MeOH'],
                                           outs=['bottom_water', 'oxygen', 'spent_catalyst', 'gas_out', 'MeOH'],
                                                water_electrolyzer=True)
    system = bst.System(path=[biorefinery_sys,
                              carbon_capture_sys])
    system.set_tolerance(mol=1e-3, rmol=1e-3, maxiter=500, subsystems=True, method='fixed point')
    return system


# def system_hydrogen_purchased(ID, **kwargs):
#     sys = create_full_MeOH_system(**kwargs)
#     @sys.flowsheet.PWC.add_specification(run=True)
#     def update_water_streams():
#         u, s = sys.flowsheet.unit, sys.flowsheet.stream
#         u.PWC.makeup_water_streams = (u.CT.ins[1], u.BT.ins[2])
#         u.PWC.process_water_streams = (
#             s.warm_process_water_1, s.ammonia_process_water,
#             s.pretreatment_steam, s.warm_process_water_2,
#             s.saccharification_water, s.stripping_water,
#             u.S401.ins[1], u.U1301.ins[2],
#             u.CIP.ins[0], u.FWT.ins[0]
#         )
#     sys.ID = ID
#     return sys

# Need to change sitepackage code of BT to exclude electrolyzer power comsumption
# IsentropicCompressor also produces power, but production = consumption;
# so change u.power_utility.consumption to rate
# self.electricity_demand = sum([u.power_utility.rate for u in units if \
                                                # u.ID != 'R1101'])
                                                
#%%

if __name__ == '__main__':
    chems = CCU.create_MeOH_chemicals()
    bst.settings.set_thermo(chems, cache=True)
    biorefinery_sys = create_cellulosic_ethanol_system('sys_ethanol_cs')
    
    emissions_BT = biorefinery_sys.flowsheet.BT.outs[0] # Boiler emissions
    emissions_fermentation = biorefinery_sys.flowsheet.D401.outs[0] # Fermentation emissions
    
    carbon_capture_sys = create_ccu_system(ins=[emissions_BT,
                                                emissions_fermentation,
                                                'makeup_MEA',
                                                'makeup_water_1',
                                                'hydrogen',
                                                'water_stream_1',
                                                'catalyst_MeOH'],
                                           outs=['bottom_water', 'oxygen', 'spent_catalyst', 'gas_out', 'MeOH'],
                                                water_electrolyzer=True)
    system = bst.System(path=[biorefinery_sys,
                              carbon_capture_sys])
    system.set_tolerance(mol=1e-5, rmol=1e-5, maxiter=1000, method='fixed-point')

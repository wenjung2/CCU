# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 10:54:48 2025

@author: IGB

"""


import biosteam as bst
import CCU

#%%
@bst.SystemFactory(ID='sys_MeOH',
                   ins=[dict(ID='emissions', CO2=1000),
                        dict(ID='makeup_MEA', H2O=1),
                        dict(ID='makeup_water_1', H2O=1),
                        dict(ID='water_stream_1', H2O=1),
                        dict(ID='catalyst_MeOH', CaSO4=1),],
                   outs=[dict(ID='bottom_water', H2O=1),
                         dict(ID='oxygen', O2=1),
                         dict(ID='spent_catalyst', CaSO4=1),
                         dict(ID='gas_out', CO2=1),
                         dict(ID='MeOH', CH3OH=1),])
def create_MeOH_system(ins, outs, water_electrolyzer=None):
    emissions, makeup_MEA, makeup_water_1, water_stream_1, catalyst_MeOH = ins
    bottom_water, oxygen, spent_catalyst, gas_out, MeOH = outs
    
    U1301 = bst.AmineAbsorption('U1301', ins=(emissions, makeup_MEA, makeup_water_1),\
                                  outs=('vent', 'concentrated'), CO2_recovery=1)
    U1301.outs[1].phase = 'g'

    ks = [bst.IsentropicCompressor(P=3*101325), 
          bst.IsentropicCompressor(P=8.8*101325),
          bst.IsentropicCompressor(P=26.3*101325)]

    hxs = [bst.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]
    
    # First 3-stage compression of CO2
    # C1101 = bst.MultistageCompressor('C1101', ins=S1301-0, outs='', compressors=ks, hxs=hxs) # won't work (no ins connection)
    C1101 = bst.MultistageCompressor('C1101', ins=U1301-1, outs='', n_stages=3, pr=2.9)
    
    # Final compression of CO2
    C1102 = bst.IsentropicCompressor('C1102', ins=C1101-0, outs='', P=78*101325)
    
    R1101 = CCU.Electrolyzer('R1101', ins=water_stream_1, outs=('', oxygen))
    @R1101.add_specification(run=True)
    def adjust_water_flow():
        C1101.run()
        R1101.ins[0].imol['H2O'] = C1101.outs[0].imol['CO2'] * 3 # H2:CO2 = 3:1
        
    C1103 = bst.IsentropicCompressor('C1103', ins=R1101-0, outs='', P=78*101325, vle=True)
        
    M1101 = bst.Mixer('M1101', ins=(C1102-0, C1103-0), outs='', rigorous=True)

    M1102 = bst.Mixer('M1102', ins=(M1101-0, ''), outs='')
    M1102.prioritize=True

    H1101 = bst.HXprocess('H1101', ins=(M1102-0, ''), outs=('', ''),)

    R1102 = CCU.MeOH_SynthesisReactor('R1102', ins=(H1101-0, catalyst_MeOH),\
                                         outs=('product', spent_catalyst))
        
    R1102_H = bst.HXutility('R1102_H', ins=R1102-0, outs='', T=284+273.15)

    S1101 = bst.Splitter('S1101', ins=R1102_H-0, outs=(1-H1101, ''), split=0.6)

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

    S1103 = bst.Splitter('S1103', ins=S1102-0, outs=('purge_2', 'recycled'), split=0.01)

    C1104 = bst.IsentropicCompressor('C1104', ins=S1103-1, outs=1-M1102, P=78*101325)

    # For condensed water and methanol in S1102 (composed of methanol, water and residual dissolved gases)
    V1102 = bst.IsenthalpicValve('V1102', ins=S1102-1, outs='', P=10*101325)

    V1103 = bst.IsenthalpicValve('V1103', ins=V1102-0, outs='', P=1.2*101325)

    F1101 = bst.SplitFlash('F1101', ins=V1103-0, outs=('gas_F1101', 1-H1103), T=22+273.15, P=1.2*101325, split=dict(CO2=1.0,
                                                                                                           H2=1.0))

    M1104 = bst.Mixer('M1104', ins=(S1103-0, F1101-0, S1104-0), outs=gas_out)
    
    T1101 = bst.StorageTank('T1101', ins=S1104-1, outs=MeOH)


def system_MeOH():
    chems = CCU.create_MeOH_chemicals()
    bst.settings.set_thermo(chems, cache=True)
    
    MeOH_sys = create_MeOH_system(ins=['CO2',
                                       'makeup_MEA',
                                       'makeup_water_1',
                                       'water_stream_1',
                                       'catalyst_MeOH'],
                                  outs=['bottom_water', 'oxygen', 'spent_catalyst', 'gas_out', 'MeOH'],
                                       )
                                                
    return MeOH_sys
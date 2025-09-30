# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 14:10:31 2025

@author: IGB
"""
import CCU 
import biosteam as bst
from numpy.testing import assert_allclose
import pytest

def test_power_consumption_after_run_with_splitting():
    system = CCU.create_full_system()
    system.set_tolerance(mol=1e-2, rmol=1e-2, subsystems=True, method='fixed point')
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = 110
    compressor = system.flowsheet.C1102
    system.run()
    power_demand_before = compressor.power_utility.consumption
    assert power_demand_before != 0, (
        'power consumption was emptied before run'
    )
    system.run()
    power_demand_after = compressor.power_utility.consumption
    assert power_demand_after != 0, (
        'power consumption was emptied after run'
    )
    assert power_demand_after > power_demand_before, (
        'expected power demand to increase after run; it did not'
    )
    
def test_power_consumption_after_convergence_with_splitting():
    system = CCU.create_full_system()
    system.set_tolerance(mol=1e-2, rmol=1e-2, subsystems=True, method='fixed point')
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = 110
    compressor = system.flowsheet.C1102
    system.run()
    power_demand_before = compressor.power_utility.consumption
    assert power_demand_before != 0, (
        'power consumption was emptied before convergence'
    )
    system.simulate()
    power_demand_after = compressor.power_utility.consumption
    assert power_demand_after != 0, (
        'power consumption was emptied after convergence'
    )
    assert power_demand_after > power_demand_before, (
        'expected power demand to increase after convergence; it did not'
    )
    
def test_emissions_accumulation_single_run():
    system = CCU.create_full_system()
    system.set_tolerance(mol=1e-2, rmol=1e-2, subsystems=True, method='fixed point')
    system.run()
    emissions_BT = system.flowsheet.emissions
    CO2_emissions_i0 = emissions_BT.imol['CO2']
    system.run()
    CO2_emissions_i1 = emissions_BT.imol['CO2']
    assert CO2_emissions_i1 > 1.2 * CO2_emissions_i0, (
         'emissions did not increase after a loop; '
        f'before {CO2_emissions_i0}, after {CO2_emissions_i1}'
    )
   
def test_emissions_convergence_no_splitting():
    system = CCU.create_full_system()
    system.set_tolerance(mol=1e-2, rmol=1e-2, subsystems=True, method='fixed point')
    system.set_tolerance(maxiter=5)
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = float('inf')
    emissions_BT = system.flowsheet.emissions
    system.run()
    CO2_emissions_0 = emissions_BT.imol['CO2']
    with pytest.raises(RuntimeError): # Should not converge
        system.simulate()
    
    CO2_emissions_f = emissions_BT.imol['CO2']
    assert CO2_emissions_f > 100 * CO2_emissions_0, (
        'system should not converge; there was no accumulation'
    )
    
def test_emissions_convergence_with_splitting():
    system = CCU.create_full_system()
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = 110
    system.set_tolerance(mol=1e-2, rmol=1e-2, subsystems=True, method='fixed point')
    system.simulate()
    CO2_emissions = splitter.outs[0].F_mol
    maximum_emissions = splitter.baseline_biorefinery_emissions * splitter.maximum_biorefinery_percent_capture / 100
    assert CO2_emissions <= maximum_emissions + 1e-64, (
        'more CO2 emissions than allowed; splitting specification is not working properly'
    )
    assert CO2_emissions > splitter.baseline_biorefinery_emissions, (
        'system emissions did not change; CCU is not working properly'
    )

    
if __name__ == '__main__':
    pass
    test_power_consumption_after_run_with_splitting()
    test_power_consumption_after_convergence_with_splitting()
    test_emissions_accumulation_single_run()
    test_emissions_convergence_no_splitting()
    test_emissions_convergence_with_splitting()
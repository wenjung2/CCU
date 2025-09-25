# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 14:10:31 2025

@author: IGB
"""
import CCU 
import biosteam as bst
from numpy.testing import assert_allclose
import pytest

def test_emissions_accumulation_single_run():
    system = CCU.create_full_system()
    emissions_BT = system.flowsheet.emissions
    CO2_emissions_i0 = emissions_BT.imol['CO2']
    system.run()
    CO2_emissions_i1 = emissions_BT.imol['CO2']
    assert CO2_emissions_i1 > 1.2 * CO2_emissions_i0, (
        f'emissions did not increase after a loop; '
        f'before {CO2_emissions_i0}, after {CO2_emissions_i1}'
    )
   
def test_emissions_convergence_no_splitting():
    system = CCU.create_full_system()
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = float('inf')
    emissions_BT = system.flowsheet.emissions
    CO2_emissions_0 = emissions_BT.imol['CO2']
    
    with pytest.raises(RuntimeError):
        system.simulate()
    
    CO2_emissions_f = emissions_BT.imol['CO2']
    assert CO2_emissions_f > 1000 * CO2_emissions_0, (
        'system should not converge; there was no accumulation'
    )
    
def test_emissions_convergence_with_splitting():
    system = CCU.create_full_system()
    emissions_BT = system.flowsheet.emissions
    splitter = system.flowsheet.carbon_capture_splitter
    splitter.maximum_biorefinery_percent_capture = 110
    CO2_emissions_0 = emissions_BT.imol['CO2']
    system.set_tolerance(rmol=0.01, subsystems=True)
    system.simulate()
    CO2_emissions_f = emissions_BT.imol['CO2']
    assert CO2_emissions_f < 1000 * CO2_emissions_0, (
        'system did not converge; there was accumulation'
    )
    
# if __name__ == '__main__':
#     test_emissions_accumulation_single_run()
#     test_emissions_convergence_no_splitting()
#     test_emissions_convergence_with_splitting()
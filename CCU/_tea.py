# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 09:58:51 2025

@author: IGB
"""


from biorefineries.cornstover import CellulosicEthanolTEA as EtOH_TEA



__all__ = ('EtOH_TEA', 'CellulosicIncentivesTEA')

#%% Incentive TEA
class CellulosicIncentivesTEA(EtOH_TEA):
    def __init__(self, *args,
                 carbon_credit=85, # $85/tonne CO2 (45Q)
                 credit_years=12,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.carbon_credit = carbon_credit
        self.credit_years = credit_years
        
    @property
    def carbon_amount_utilized(self):
        if "MeOH" not in self.system.flowsheet.stream:
            return 0.0 
        else:
            methanol = self.system.flowsheet.MeOH
            hours = self.operating_days * 24
            CO2_utilized = methanol.get_atomic_flow('C') * 32.04 * hours/1000 # in ton/year
            return CO2_utilized

    @property
    def annual_credit(self):
        return self.carbon_amount_utilized * self.carbon_credit
    
    def _fill_tax_and_incentives(self, incentives, taxable_cashflow, nontaxable_cashflow, tax, depreciation):
        super()._fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax, depreciation)
        if self.carbon_amount_utilized is None:
            return
        annual_credit = (self.carbon_amount_utilized * self.carbon_credit)
        
        for year in range(len(incentives)):
            if year < 3:
                continue
            
            project_year = year - 3 + 1
            if project_year <= self.credit_years:
                incentives[year]  += annual_credit

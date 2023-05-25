from AuxiliaryComponent import AuxiliaryComponent
from Constants import *


class BiomassGasifier(AuxiliaryComponent):
    def __init__(self, id, replacement_year,
                 biomass_flow_rate, eff,
                 en_el_abs_rate, heat_developed_rate, biomass, syngas, temperature_syngas,
                 pressure_syngas, replacement_cost, tech='converter',heat_recovery_system=None):
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost)

        """
        :param tech:'storage', 'pump','compressor', 'converter'
        :param id: identification code
        :param replacement_year [year]
        :param biomass_flow_rate:biomass nominal flow rate [kg/h]
        :param eff: cold gas efficiency: energy of purifed syngas/ energy of biomass
        :param en_el_abs_rate: electric energy absorbed per kg of biomass [kWh/kg]
        :param heat_developed_rate: heat power recovered in the syngas purification phase per kg of biomass [kWh/kg]
        :param biomass: type of input biomass (wood pellets,wood chips,other)
        :param syngas: typer of output syngas 
        :param temperature_syngas: purification process downstream temperature [°C]
        :param pressure_syngas: purification process downstream pressure [bar]
        """

        self.biomass_flow_rate = biomass_flow_rate
        self.eff = eff
        self.en_el_abs_rate = en_el_abs_rate
        self.q_recovered_rate = heat_developed_rate
        self.biomass = biomass
        self.syngas = syngas
        self.temperature_syngas = temperature_syngas
        self.pressure_syngas = pressure_syngas
        self.heat_recovery_system=heat_recovery_system

    def compute_output(self, time):
        """
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :return: syngas_flow_rate: [kg/time]
        :return: syngas_flow_rate_1: [Nm3/time]
        :return: power_syngas: [kW]
        :return: power_el_abs: [kW]
        :return: power_heat_developed: [kW]
        """

        syngas_flow_rate = self.eff * self.biomass_flow_rate * time * self.biomass.lhv / self.syngas.lhv

        energy_syngas = syngas_flow_rate * self.syngas.lhv
        power_syngas = energy_syngas / (conv_hour_sec * time)

        en_el_abs = self.en_el_abs_rate * self.biomass_flow_rate * time
        heat_developed = self.q_recovered_rate * self.biomass_flow_rate * time

        biomass_cons = self.biomass_flow_rate * time
        power_el_abs = en_el_abs / time
        power_heat_developed = heat_developed / time
        #syngas_flow_rate_1 = (syngas_flow_rate * r_gas * (t_ref + conv_celsius_kelvin)) / (
         #       self.syngas.molecular_weight/conv_g_kg * p_ref * conv_bar_pascal)
        syngas_flow_rate_1=syngas_flow_rate/self.syngas.rho
        if self.heat_recovery_system:
            power_heat_recovered = power_heat_developed*self.heat_recovery_system.eff
        else:
            power_heat_recovered=0



        return biomass_cons, syngas_flow_rate, syngas_flow_rate_1, power_syngas, power_el_abs, power_heat_developed,power_heat_recovered

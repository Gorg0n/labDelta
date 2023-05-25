from AuxiliaryComponent import AuxiliaryComponent
from Constants import *
import numpy as np


class GasStorage(AuxiliaryComponent):
    def __init__(self, id, capacity,
                 pressure_max, replacement_year, replacement_cost, tech='GasStorage'):
        """

        :param id:
        :param capacity: tank volume [m3]
        :param pressure_max: maximum admissible pressure [bar]
        :param replacement_year:
        :param replacement_cost: initial cost [€]
        :param tech:
        """
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost)

        self.capacity = capacity
        self.pressure_max = pressure_max

    def compute_output(self, gas_flow_rate_input, gas_flow_rate_output, initial_pressure_level, temperature_storage):
        """

        :param gas_flow_rate_input: [Nm3/time]
        :param gas_flow_rate_output: [Nm3/time]
        :param initial_pressure_level: initial pressure/maximum allowable pressure ∈[0,1]
        :param temperature_storage: temperature of gas[°C]
        :return:
        pressure: pressure in the storage [bar]
        pressure_level: pressure/maximum allowable pressure ∈[0,1]
        v: volume stored in tank [Nm3]
        v_surplus: surplus volume of gas [Nm3/time]
        v_unmet: deficit volume of gas [Nm3/time]

        """
        pressure = np.zeros((len(temperature_storage),))
        pressure_level = np.zeros((len(temperature_storage),))
        v = np.zeros((len(temperature_storage),))
        v_surplus = np.zeros((len(temperature_storage),))
        v_unmet = np.zeros((len(temperature_storage),))

        rho_ref = p_ref * conv_bar_pascal / (r_gas * (t_ref + conv_celsius_kelvin))
        p_initial = initial_pressure_level * self.pressure_max * conv_bar_pascal
        mol = p_initial * self.capacity / (r_gas * (temperature_storage[0] + conv_celsius_kelvin))
        for i in range(len(temperature_storage)):
            temperature_storage_i = temperature_storage[i] + conv_celsius_kelvin
            mol_max = self.pressure_max * conv_bar_pascal * self.capacity / (r_gas * temperature_storage_i)
            mol_storable = mol_max - mol
            mol_available = mol
            mol_input = gas_flow_rate_input[i] * rho_ref
            mol_output = gas_flow_rate_output[i] * rho_ref

            mol_balance = mol_input - mol_output

            if mol_balance >= 0:
                if mol_balance >= mol_storable:
                    mol = mol_max
                    mol_surplus = mol_balance - mol_storable
                    mol_unmet = 0
                else:
                    mol += mol_balance
                    mol_surplus = 0
                    mol_unmet = 0
            else:
                if abs(mol_balance) >= mol_available:
                    mol = 0
                    mol_surplus = 0
                    mol_unmet = abs(mol_balance) - mol_available
                else:
                    mol -= abs(mol_balance)
                    mol_surplus = 0
                    mol_unmet = 0

            pressure_i = (mol * r_gas * temperature_storage_i / self.capacity) / conv_bar_pascal
            pressure_level_i = pressure_i / self.pressure_max

            pressure[i] = pressure_i
            pressure_level[i] = pressure_level_i
            v[i] = mol / rho_ref
            v_surplus[i] = mol_surplus / rho_ref
            v_unmet[i] = mol_unmet / rho_ref


        return pressure, pressure_level, v, v_surplus, v_unmet

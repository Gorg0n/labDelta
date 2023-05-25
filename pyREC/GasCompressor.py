from AuxiliaryComponent import AuxiliaryComponent
from Constants import *
import numpy as np


class GasCompressor(AuxiliaryComponent):
    def __init__(self, id, fuel, replacement_year, replacement_cost, tech='compressor', ):
        """

        :param id:
        :param tech:
        :param replacement_year:
        :param fuel: gas fuel
        """
        super().__init__(id=id, tech=tech, replacement_year=replacement_year, replacement_cost=replacement_cost)

        self.fuel = fuel

    def compute_output(self, time, temperature_gas_input, pressure_gas_input, gas_flow_rate_input, pressure_gas_output):
        """
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :param temperature_gas_input: [°C]
        :param pressure_gas_input: [bar]
        :param gas_flow_rate_input: [Nm3/time]
        :param pressure_gas_output: [bar]
        :return:
        power_el_abs: [W]
        q_power_developed: [W]
        temperature_gas_output: [°C]
        """
        n_stages=5
        n_parallel=1
        power_el_abs = np.zeros((len(temperature_gas_input),))
        q_power_developed = np.zeros((len(temperature_gas_input),))
        temperature_gas_output = np.zeros((len(temperature_gas_input),))

        rho_ref = p_ref * conv_bar_pascal / (r_gas * (t_ref + conv_celsius_kelvin))
        c_v = self.fuel.c_p - r_gas
        polytropic_coeff = self.fuel.c_p / c_v
        for i in range(len(temperature_gas_input)):
            temperature_gas_input_i = temperature_gas_input[i] + conv_celsius_kelvin
            pressure_ratio = (pressure_gas_output[i] / pressure_gas_input[i])**(1/n_stages)
            p_stages=np.zeros([n_stages+1,1])
            p_stages[0]=pressure_gas_input[i]
            p_stages[n_stages]=pressure_gas_output[i]
            for stage in range(1,n_stages):
                p_stages[stage]=pressure_ratio**(stage)*p_stages[0]
            w_stage=(polytropic_coeff*r_gas*temperature_gas_input_i)/(polytropic_coeff-1)*(1-(p_stages[1]/p_stages[0])**((polytropic_coeff-1)/polytropic_coeff))
            w_polytropic=w_stage*n_stages
            mol_input = gas_flow_rate_input[i] * rho_ref
            power_el_abs[i] = n_parallel*abs(w_polytropic) * mol_input / (conv_hour_sec * time)
            temperature_gas_output_i = (
                        temperature_gas_input_i * (pressure_ratio) ** ((polytropic_coeff - 1) / polytropic_coeff))
            q_power_developed[i] = mol_input / (conv_hour_sec * time) * self.fuel.c_p * (
                        temperature_gas_output_i - temperature_gas_input_i)*n_stages*n_parallel
            temperature_gas_output[i] = temperature_gas_output_i - conv_celsius_kelvin



            # temperature_gas_input_i = temperature_gas_input[i] + conv_celsius_kelvin
            # pressure_ratio = pressure_gas_output[i] / pressure_gas_input[i]
            # w_polytropic = ((polytropic_coeff * r_gas * temperature_gas_input_i) / (polytropic_coeff - 1)) * (
            #             1 - (pressure_ratio) ** ((polytropic_coeff - 1) / polytropic_coeff))
            # mol_input = gas_flow_rate_input[i] * rho_ref
            # power_el_abs[i] = abs(w_polytropic) * mol_input / (conv_hour_sec * time)
            # temperature_gas_output_i = (
            #             temperature_gas_input_i * (pressure_ratio) ** ((polytropic_coeff - 1) / polytropic_coeff))
            # q_power_developed[i] = mol_input / (conv_hour_sec * time) * self.fuel.c_p * (
            #             temperature_gas_output_i - temperature_gas_input_i)
            # temperature_gas_output[i] = temperature_gas_output_i - conv_celsius_kelvin

        return power_el_abs, q_power_developed, temperature_gas_output

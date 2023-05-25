from ProductionSystem_deleate import ProductionSystem
from Fuel import Fuel
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.interpolate import interpn


class GasTurbine(ProductionSystem):
    def __init__(self, id, cp=None, tech='GT', user=None, status='new', power=1,
                 decay=0.006,
                 life_time=20, cost_kW=1500, oem_cost_kW=40, inc=None, dur_inc=None,
                 storage=None, mode=None, datasheet_curves_ref=None, datasheet_curves=None, fuel=None):
        super().__init__(id=id, cp=cp, tech=tech, user=user, status=status, power=power,
                         decay=decay,
                         life_time=life_time, cost_kW=cost_kW, oem_cost_kW=oem_cost_kW, inc=inc,
                         dur_inc=dur_inc,
                         storage=storage)

        """
        :param id: identification code
        :param cp:
        :param tech: 'gas turbine'
        :param user: list of user physically connected
        :param status: 'new' or 'old'
        :param power: rated power [kW]
        :param decay:decay of production plants (for example 0.006)
        :param life_time:useful life of the plant [year]
        :param cost_kW: initial cost [€/kW]
        :param oem_cost_kW:maintenance costs [€/year/kW]
        :param inc: type of incentive for example bonus50%
        :param dur_inc: duration incentive [year]
        :param storage:
        :param mode:1 datasheet_curves, 2 datasheet_curves_ref 
        :param  datasheet_curves : multi index dataframe including variables dependence on temperature and ratio_power, from file.xlsx pd.read_excel('file_path/file_name', header=[0, 1]
        :param datasheet_curves_ref: multi index dataframe including variables dependence on temperature and ratio_power, from file.xlsx pd.read_excel('file_path/file_name', header=[0, 1]
        :param fuel: LPG (liquefied gas) , propane (C3H8), methane (CH4), natural gas, hydrogen (H2), syngas, biogas,
        
        """
        self.mode = mode
        self.datasheet_curves = datasheet_curves
        self.datasheet_curves_ref = datasheet_curves_ref
        self.fuel = fuel

        if self.mode == 1:
            self.df = self.datasheet_curves
        else:
            self.df = self.datasheet_curves_ref

    def power_temperature_curve(self):
        """

        :return: f interpolator function power_max=f(T_amb)
        """
        df0 = self.df['Capacity [kW]']
        power_max_0 = df0.iloc[len(df0) - 1, :]
        t0 = np.array(list(df0.columns.values))
        f_p_max = interp1d(t0, power_max_0, kind='cubic')
        return f_p_max

    def datasheet_interpolation(self, power_ratio, t_amb, ref_values=None):
        """
        :param: power_ratio: Power/Max Power
        :param: t_amb:ambient temperature [°C]
        :param: list of reference_values at T=20° and ratio_power=1  [Capacity [kW],Heat Rate [kJ/kWh],Exhaust Temperature [C],Exhaust Heat [GJ/hr],Air Inlet Flow [kg/hr],Exhaust Flow Rate [kg/hr]]
        :return:
        """

        power_ratio_0 = np.array(self.df.iloc[:, 0])

        df0 = self.df['Capacity [kW]']
        len_ref = len(list(df0.columns.values))
        t0 = np.array(list(df0.columns.values))

        header = []
        for i, k in enumerate(self.df.columns.tolist()[1:]):
            e = k[0]
            header.append(e)
        header = list(set(header))
        variables = np.zeros(len(header))

        for k in range(len(header)):
            df2 = self.df.iloc[:, 1:]
            df2 = df2.iloc[:, k * len_ref:(k + 1) * len_ref]
            values = np.array(df2)
            variables[k] = interpn((power_ratio_0, t0), values, np.array([power_ratio, t_amb]))

        if self.mode == 1:
            variables = variables
        else:
            variables = variables * np.array(ref_values)

        return variables

    def compute_output(self, t_amb, d, time, ref_values=None):
        """

        :param t_amb:ambient temperature [°C]
        :param d:electricity power demand [kW]
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :param: list of reference values at T=20° and ratio_power=1: [Capacity [kW],Heat Rate [kJ/kWh],Exhaust Temperature [C],Exhaust Heat [GJ/hr],Air Inlet Flow [kg/hr],Exhaust Flow Rate [kg/hr]]
        :return:
        """

        f_power_max = self.power_temperature_curve()
        if self.mode == 1:
            power_max = f_power_max(t_amb)
        else:
            power_max_normalized = f_power_max(t_amb)
            power_max = power_max_normalized * ref_values[0]

        if d >= power_max:
            p = power_max
            surplus = 0
            p_unmet = d - p
            power_ratio = 1


        elif d < power_max and d != 0:
            p = d
            surplus = 0
            p_unmet = 0
            power_ratio = d / power_max

        else:
            p = 0
            surplus = 0
            p_unmet = 0
            power_ratio = 0

        if power_ratio < 0.1 and power_ratio != 0:
            power_ratio = 0.1

        variables = self.datasheet_interpolation(power_ratio=power_ratio, t_amb=t_amb, ref_values=ref_values)
        q_abs_rate = variables[1]
        q_wasted_rate = variables[2]
        t_out = variables[3]
        air_in_rate = variables[4]
        exhaust_flow_rate = variables[5]

        q_abs = q_abs_rate * p * time
        q_wasted = q_wasted_rate * time
        air_in = air_in_rate * time
        exhaust_flow = exhaust_flow_rate * time
        fuel_cons = exhaust_flow - air_in

        fuel_cons_calc = q_abs / self.fuel.lhv
        if self.fuel.state == 'liquid':
            tCO2 = fuel_cons_calc * self.fuel.rho_liq / 1000 * self.fuel.tCO2_emission
        elif self.fuel.state == 'solid':
            tCO2 = fuel_cons_calc / 1000 * self.fuel.tCO2_emission
        elif self.fuel.state == 'gas':
            tCO2 = fuel_cons_calc * self.fuel.rho_gas / 1000 * self.fuel.tCO2_emission
        else:
            tCO2 = None
            print("lhv and rho must be consistent with unit_fuel_cons")

        return power_ratio, p, surplus, p_unmet, q_abs_rate, q_wasted_rate, t_out, air_in_rate, exhaust_flow_rate, fuel_cons, fuel_cons_calc, tCO2


datasheet_curves_ref = pd.read_excel('DataSheetGasTurbine_reference.xlsx', header=[0, 1])

fuel=Fuel(fuel='LGP')
gas_turbine1 = GasTurbine(id='ciao', cp=None, tech='GT', user=None, status='new', power=1,
                          decay=0.006,
                          life_time=20, cost_kW=1500, oem_cost_kW=40, inc=None, dur_inc=None,
                          storage=None, mode=1,
                          datasheet_curves_ref=pd.read_excel('DataSheetGasTurbine_reference.xlsx', header=[0, 1]),
                          datasheet_curves=pd.read_excel('DataSheetGasTurbine.xlsx', header=[0, 1]), fuel=fuel)

output = gas_turbine1.compute_output(t_amb=20, d=10, time=0.25, ref_values=[3387, 12995, 440.1, 29.931, 64567, 65445])

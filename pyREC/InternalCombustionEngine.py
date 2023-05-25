from ProductionSystem import ProductionSystem
from AuxiliaryComponent import AuxiliaryComponent
import math
import numpy as np
from scipy.interpolate import interp1d
from Constants import *


class InternalCombustionEngine(ProductionSystem, AuxiliaryComponent):
    def __init__(self, id, carrier, power, category, datasheet_curve, power_max, power_min, n_min,
                 n_max,
                 fuel, cogeneration,
                 p_con, cost_kW, inc_kW, oem_cost_kW, oem_cost_kWh, inc_kWh, inc, dur_inc_kWh,
                 dur_inc, user=None, pod=None, ma=None, ts=None, gse_mode='rid',
                 tech='ICE', status='new', decay=0.006,
                 life_time=20, aux_components=None,
                 replacement_year=None, replacement_cost=None):
        ProductionSystem.__init__(self, id=id, carrier=carrier, tech=tech, power=power, category=category,
                                  status=status,
                                  user=user, decay=decay,
                                  life_time=life_time, pod=pod, ma=ma, ts=ts, gse_mode=gse_mode,
                                  p_con=p_con, cost_kW=cost_kW, inc_kW=inc_kW, inc_kWh=inc_kWh, oem_cost_kW=oem_cost_kW,
                                  oem_cost_kWh=oem_cost_kWh,
                                  inc=inc, dur_inc=dur_inc, dur_inc_kWh=dur_inc_kWh,
                                  aux_components=aux_components, cogeneration=cogeneration)

        AuxiliaryComponent.__init__(self, id=id, tech=None,
                                    replacement_year=replacement_year, replacement_cost=replacement_cost)

        """
        :param id: identification code
        :param cp:
        :param tech:InternalCombustionEngine (ICE)
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
        :param replace: (True or False)
        :param replacement_year
        :param mode: 1 based on demand (ICE is flexible) 2 based on fuel flow rate
        :param datasheet_curve: fuel consumption [unit fuel cons*/h] *l/h  or Nm3/h vs Power/Rated Power
        :param power_max: maximum output power [kW]
        :param power_min: minimum output power [kW]
        :param n_min: minimum number of ICE in operation
        :param n_max: maximum number of ICE in operation 
        :param fuel: LPG (liquefied gas) , propane (C3H8), methane (CH4), natural gas, hydrogen (H2), syngas, biogas
     
        """
        if self.category == 'flex':
            self.mode = 1
        else:
            self.mode = 2
        self.datasheet_curve = datasheet_curve
        self.power_max = power_max
        self.power_min = power_min
        self.n_min = n_min
        self.n_max = n_max
        self.fuel = fuel
        self.cogeneration = cogeneration

    def fuel_consumption_curve(self):
        """
        :return: fuel consumption curve
        """
        x = np.array(self.datasheet_curve['power_ratio'])
        y = np.array(self.datasheet_curve['fuel_cons'])
        f = interp1d(x, y, kind='cubic')
        f1 = interp1d(y, x, kind='cubic')
        if self.mode == 1:
            return f
        else:
            self.fuel_flow_rate_min = f(x[0])
            self.fuel_flow_rate_max = f(x[len(x) - 1])
            return f1

    def controller(self, time=None, fuel_flow_rate=None, d=None):
        """
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis input in mode 2
        :param d: electricity power demand [kW] input in mode 1
        :param fuel_flow_rate: [unit fuel/unit time] input in mode 2
        :return: n: Total number of ICE required to meet load.
        :return: p_set: power set point [kW]
        :return: p: [kW]
        :return surplus: [kW]
        :return p_unmet:[kW]
        """

        if self.mode == 1:
            n_ideal = math.ceil(d / self.power)

            if self.n_max == self.n_min:
                n = self.n_max
                p_set_ideal = d / n
                if p_set_ideal >= self.power_max:
                    p_set = self.power_max
                elif p_set_ideal <= self.power_min:
                    p_set = self.power_min
                else:
                    p_set = p_set_ideal
            else:
                if n_ideal >= self.n_max:
                    n = self.n_max
                    p_set_ideal = d / self.n_max
                    if p_set_ideal >= self.power_max:
                        p_set = self.power_max
                    else:
                        p_set = p_set_ideal
                elif n_ideal <= self.n_min:
                    n = self.n_min
                    if n == 0:
                        p_set_ideal = 0
                    else:
                        p_set_ideal = d / n
                    if p_set_ideal <= self.power_min:
                        p_set = self.power_min
                    else:
                        p_set = p_set_ideal
                else:
                    n = n_ideal
                    p_set_ideal = d / n_ideal
                    p_set = p_set_ideal

            p = p_set * n
            if p >= d:
                surplus = p - d
                unmet = 0
                self_cons = d
            else:
                surplus = 0
                unmet = d - p
                self_cons=p

            return n, p_set, p, surplus, unmet,self_cons

        else:
            f = self.fuel_consumption_curve()
            fuel_flow_rate_min = self.fuel_flow_rate_min
            fuel_flow_rate_max = self.fuel_flow_rate_max
            fuel_flow_rate_h = fuel_flow_rate / time
            n_ideal = math.ceil(fuel_flow_rate_h / fuel_flow_rate_max)

            if self.n_max == self.n_min:
                n = self.n_max
                fuel_flow_rate_set_h = fuel_flow_rate_h / n
                if fuel_flow_rate_set_h > fuel_flow_rate_max:
                    fuel_flow_rate_set_h = fuel_flow_rate_max
                    fuel_flow_surplus_h = fuel_flow_rate_h - n * fuel_flow_rate_set_h
                    fuel_flow_deficit_h = 0

                elif fuel_flow_rate_set_h < fuel_flow_rate_min:
                    print('min numero di motori accessi non può essere garantito per mancanza combustibile')
                    while n > 0 and fuel_flow_rate_set_h < fuel_flow_rate_min:
                        n -= 1
                        if n != 0:
                            fuel_flow_rate_set_h = fuel_flow_rate_h / n
                        else:
                            fuel_flow_rate_set_h = 0

                    if fuel_flow_rate_set_h > fuel_flow_rate_max:
                        fuel_flow_rate_set_h = fuel_flow_rate_max
                    else:
                        fuel_flow_rate_set_h = fuel_flow_rate_set_h

                    fuel_flow_deficit_h = self.n_min * fuel_flow_rate_min - fuel_flow_rate_h
                    fuel_flow_surplus_h = 0
                else:
                    fuel_flow_rate_set_h = fuel_flow_rate_set_h
                    fuel_flow_deficit_h = 0
                    fuel_flow_surplus_h = 0

            else:
                if n_ideal > self.n_max:
                    n = self.n_max
                    fuel_flow_rate_set_h = fuel_flow_rate_max
                    fuel_flow_surplus_h = fuel_flow_rate_h - n * fuel_flow_rate_set_h
                    fuel_flow_deficit_h = 0

                else:
                    if n_ideal < self.n_min:
                        n = self.n_min
                    else:
                        n = n_ideal

                    fuel_flow_rate_set_h = fuel_flow_rate_h / n
                    if fuel_flow_rate_set_h > fuel_flow_rate_max:
                        fuel_flow_rate_set_h = fuel_flow_rate_max
                        fuel_flow_surplus_h = fuel_flow_rate_h - n * fuel_flow_rate_set_h
                        fuel_flow_deficit_h = 0

                    elif fuel_flow_rate_set_h < fuel_flow_rate_min:
                        while n > 0 and fuel_flow_rate_set_h < fuel_flow_rate_min:
                            n -= 1
                            if n != 0:
                                fuel_flow_rate_set_h = fuel_flow_rate_h / n
                            else:
                                fuel_flow_rate_set_h = 0

                        if n < self.n_min and fuel_flow_rate_set_h > fuel_flow_rate_max:
                            print('min numero di motori accessi non può essere garantito per mancanza combustibile')
                            fuel_flow_rate_set_h = fuel_flow_rate_max
                            fuel_flow_deficit_h = self.n_min * fuel_flow_rate_min - fuel_flow_rate_h
                            fuel_flow_surplus_h = 0
                        elif n < self.n_min and fuel_flow_rate_set_h < fuel_flow_rate_max:
                            print('min numero di motori accessi non può essere garantito per mancanza combustibile')
                            fuel_flow_rate_set_h = fuel_flow_rate_set_h
                            fuel_flow_deficit_h = self.n_min * self.fuel_flow_rate_min - fuel_flow_rate_h
                            fuel_flow_surplus_h = 0

                        elif n > self.n_min and fuel_flow_rate_set_h > fuel_flow_rate_max:
                            fuel_flow_rate_set_h = self.fuel_flow_rate_max
                            fuel_flow_surplus_h = fuel_flow_rate_h - fuel_flow_rate_set_h * n
                            fuel_flow_deficit_h = 0
                        else:
                            fuel_flow_rate_set_h = fuel_flow_rate_set_h
                            fuel_flow_surplus_h = fuel_flow_rate_h - fuel_flow_rate_set_h * n
                            fuel_flow_deficit_h = 0


                    else:
                        fuel_flow_rate_set_h = fuel_flow_rate_set_h
                        fuel_flow_surplus_h = 0
                        fuel_flow_deficit_h = 0

            if fuel_flow_rate_set_h < fuel_flow_rate_min:
                ratio_power = 0
            else:
                ratio_power = f(fuel_flow_rate_set_h)

            p_set = ratio_power * self.power

            p = p_set * n

            fuel_flow_rate_set = fuel_flow_rate_set_h * time
            fuel_flow_surplus = fuel_flow_surplus_h * time
            fuel_flow_deficit = fuel_flow_deficit_h * time

            return n, p_set, p, fuel_flow_rate_set, fuel_flow_surplus, fuel_flow_deficit

    def compute_output(self, time, d=None, fuel_flow_rate=None):
        """
        :param time: 1 if hourly analysis [h], 0.25 if quarterly analysis [1/4 h]
        :param d: electricity power demand [kW], (or heat demand [kW] if carrier='heat')
        :param fuel_flow_rate: [unit fuel/time]
        :return: n: Total number of ICE required to meet load.
        :return: p_set: power set point [kW]
        :return: p_tot: [kW]
        :return surplus: [kW]
        :return unmet:[kW]
        :return power_heat_developed:[kW]
        :return fuel consumption rate: [unit fuel cons/time]
        :return:power_heat_developed: [kW]
        :return:eta_fuel: [kWh/(unit fuel cons)]
        :return:eta_el: Electrical efficiency
        :return:fuel_cons: liquid fuel consumption [unit fuel cons/time]
        :return:tC02_tot: tC02 emission [ton/time]

        """
        if d is not None:
            len_ref = len(d)
        else:
            len_ref = len(fuel_flow_rate)

        n = np.zeros((len_ref,))
        p_set = np.zeros((len_ref,))
        p = np.zeros((len_ref,))
        surplus = np.zeros((len_ref,))
        unmet = np.zeros((len_ref,))
        self_cons = np.zeros((len_ref,))
        power_heat_developed = np.zeros((len_ref,))
        power_heat_recovered = np.zeros((len_ref,))
        eta_fuel = np.zeros((len_ref,))
        eta_el = np.zeros((len_ref,))
        fuel_cons = np.zeros((len_ref,))
        tCO2 = np.zeros((len_ref,))
        fuel_flow_surplus = np.zeros((len_ref,))
        fuel_flow_deficit = np.zeros((len_ref,))
        surplus_heat = np.zeros((len_ref,))
        unmet_heat = np.zeros((len_ref,))
        self_cons_heat = np.zeros((len_ref,))

        for i in range(len_ref):
            if self.mode == 1:
                n_i, p_set_i, p_i, surplus_i, unmet_i,self_cons_i = self.controller(d=d[i])
                if self.carrier == 'heat':
                    eta_el_ref = 0.35
                    n_i, p_set_i, p_i, surplus_i, unmet_i,self_cons_i = self.controller(d=d[i] * eta_el_ref / (1 - eta_el_ref)*self.aux_components[0].eff)
                x = p_set_i / self.power
                if x == 0:
                    y = 0

                else:
                    f = self.fuel_consumption_curve()
                    y = f(x) * time
                fuel_flow_surplus_i = 0
                fuel_flow_deficit_i = 0
            else:
                n_i, p_set_i, p_i, fuel_flow_rate_set_i, fuel_flow_surplus_i, fuel_flow_deficit_i = self.controller(
                    fuel_flow_rate=fuel_flow_rate[i], time=time)
                y = fuel_flow_rate_set_i
                surplus_i = None
                unmet_i = None
                self_cons_i=None

            if p_set_i == 0:
                eta_fuel_i = 0
                eta_el_i = 0
                tCO2_i = 0
                power_heat_developed_i = 0
            else:
                eta_fuel_i = p_set_i * time / y  # [kWh/unit fuel]
                eta_el_i = p_set_i * conv_hour_sec * time / (self.fuel.rho * self.fuel.lhv * y)
                #eta_el_i=0.35 #attenzione dipende fortemente dai consumii
                print(eta_el_i)
                tCO2_i = y * self.fuel.rho / 1000 * self.fuel.tCO2_emission
                power_heat_developed_i = p_set_i * (1 - eta_el_i) / eta_el_i

            power_heat_developed_i = n_i * power_heat_developed_i
            fuel_cons_i = n_i * y
            tCO2_i = n_i * tCO2_i

            if self.carrier == 'heat':
                p_heat = power_heat_developed_i*self.aux_components[0].eff
                d_heat = d[i]
                if p_heat >= d_heat:
                    surplus_heat_i = p_heat - d_heat
                    deficit_heat_i = 0
                    self_cons_heat_i=d_heat
                else:
                    surplus_heat_i = 0
                    deficit_heat_i = d_heat-p_heat
                    self_cons_heat_i = p_heat

                surplus_heat[i] = surplus_heat_i
                unmet_heat[i] = deficit_heat_i
                self_cons_heat[i]=self_cons_heat_i

            n[i] = n_i
            p_set[i] = p_set_i
            p[i] = p_i
            surplus[i] = surplus_i
            unmet[i] = unmet_i
            self_cons[i]=self_cons_i
            power_heat_developed[i] = power_heat_developed_i
            eta_fuel[i] = eta_fuel_i
            eta_el[i] = eta_el_i
            fuel_cons[i] = fuel_cons_i
            tCO2[i] = tCO2_i
            fuel_flow_surplus[i] = fuel_flow_surplus_i
            fuel_flow_deficit[i] = fuel_flow_deficit_i
            if self.cogeneration==True:
                power_heat_recovered[i] = power_heat_developed_i*self.aux_components[0].eff


        self.en_perf_evolution['electricity'] = {}
        self.en_perf_evolution['heat'] = {}
        self.en_perf_evolution['electricity']['prod'] = p
        self.en_perf_evolution['electricity']['surplus'] = surplus
        self.en_perf_evolution['electricity']['unmet'] = unmet
        self.en_perf_evolution['electricity']['self_cons'] =self_cons
        self.en_perf_evolution['heat']['prod'] = power_heat_recovered

        if self.carrier == 'heat':
            self.en_perf_evolution['heat']['surplus'] = surplus_heat
            self.en_perf_evolution['heat']['unmet'] = unmet_heat
            self.en_perf_evolution['heat']['self_cons'] = self_cons_heat

        self.cost_resource = sum(fuel_cons) * self.fuel.cost
        self.fuel_cons = fuel_cons

        self.output=n, p_set, p, surplus, unmet, power_heat_developed,power_heat_recovered, eta_fuel, eta_el, fuel_cons, tCO2, fuel_flow_surplus, fuel_flow_deficit
        return n, p_set, p, surplus, unmet, power_heat_developed,power_heat_recovered, eta_fuel, eta_el, fuel_cons, tCO2, fuel_flow_surplus, fuel_flow_deficit


    def check_eff_calculate_cb(self, time, value_cb,self_cons_heat=None):

        fuel_cons = self.fuel_cons  # Nm3
        power_fuel = fuel_cons * self.fuel.rho * self.fuel.lhv / (3600 * time)  # The produced fuel power (kW)

        prod_el = self.en_perf_evolution['electricity']['prod']
        if self.carrier=='heat':
            self_cons_heat = self.en_perf_evolution['heat']['self_cons'] * self.aux_components[0].eff
        else:
            self_cons_heat=self_cons_heat


        h = 0
        for element in prod_el:
            if element > 0:
                h += 1

        eff_heat = self_cons_heat / power_fuel
        print(eff_heat,'eff_heat')
        eff_el = prod_el / power_fuel
        print(eff_el,'eff_el')
        eff_heat_ref_list = [0.89, 0.89, 0.9, 0.9, 0.89, 0.7]  # 0.90- 0.70 methane and syngas
        eff_el_ref_list = [0.442, 0.442, 0.525, 0.525, 0.442, 0.42]  # 0.525-0.42
        if self.fuel.fuel == 'LGP':
            eff_heat_ref = eff_heat_ref_list[0]
            eff_el_ref = eff_el_ref_list[0]
        elif self.fuel.fuel == 'diesel':
            eff_heat_ref = eff_heat_ref_list[1]
            eff_el_ref = eff_el_ref_list[1]

        elif self.fuel.fuel == 'methane':
            eff_heat_ref = eff_heat_ref_list[2]
            eff_el_ref = eff_el_ref_list[2]

        elif self.fuel.fuel == 'natural gas':
            eff_heat_ref = eff_heat_ref_list[3]
            eff_el_ref = eff_el_ref_list[3]

        elif self.fuel.fuel == 'hydrogen':
            eff_heat_ref = eff_heat_ref_list[4]
            eff_el_ref = eff_el_ref_list[4]

        elif self.fuel.fuel == 'syngas':
            eff_heat_ref = eff_heat_ref_list[5]
            eff_el_ref = eff_el_ref_list[5]

        pes = 1 - 1 / (eff_heat / eff_heat_ref + eff_el / eff_el_ref)

        pes_car = np.zeros((len(pes),))
        for dt, element in enumerate(pes):
            if element > 0:
                pes_car[dt] = 1
            else:
                pes_car[dt] = 0

        eff_el_ref_1 = 0.46

        risp = (prod_el / eff_el_ref_1 + self_cons_heat / eff_heat_ref - power_fuel) * time / 1000  # MWh
        for dt,element in enumerate(risp):
            if element<0:
                risp[dt]=0
            else:
                risp[dt]=element

        # # To obtain the vale of K, we need to have this ratio:
        # # q = The amount of energy produced in CAR criteria / number of working hours of CHP

        q = sum(prod_el * pes_car) / h  # kW

        if q <= 1:
            k = 1.4 * 1
        elif 1 < q <= 10:
            k = 1.4 * 1 + (q - 1) * 1.3
        elif 10 < q <= 80:
            k = 1.4 * 1 + 9 * 1.3 + (q - 10) * 1.2
        elif 80 < q <= 100:
            k = 1.4 * 1 + 9 * 1.3 + 70 * 1.2 + (q - 80) * 1.1
        elif q > 100:
            k = 1 * q
        k = k / q
        cb = sum(risp) * 0.086 * k  # MWh
        cbt = cb * value_cb  # noticing that every Certificato Bianco values between 250 and 260 Euros/MWh

        return pes, pes_car, h, cbt

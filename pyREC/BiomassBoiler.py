from ProductionSystem import ProductionSystem
import numpy as np
from scipy.interpolate import interp1d
from Constants import *


class BiomassBoiler(ProductionSystem):
    def __init__(self, id, power, cost_kW, inc_kW, oem_cost_kW, inc_kWh, inc, dur_inc_kWh,
                 dur_inc, power_el_absorbed, eff,
                 biomass_flow_rate_max, biomass_flow_rate_min, biomass, category, carrier='heat', tech='biomass_boiler',
                 status='new', p_con=None, user=None, decay=0.006,
                 life_time=20,
                 aux_components=None, pm_emmissiosn_rate=None,
                 co_emissions_rate=None, cogeneration=False):
        super().__init__(id=id, carrier=carrier, tech=tech, power=power, status=status, user=user, decay=decay,
                         life_time=life_time,
                         p_con=p_con, cost_kW=cost_kW, inc_kW=inc_kW, inc_kWh=inc_kWh, oem_cost_kW=oem_cost_kW, inc=inc,
                         dur_inc=dur_inc, dur_inc_kWh=dur_inc_kWh,
                         aux_components=aux_components, category=category, cogeneration=cogeneration)

        """
        :param id: identification code
        :param carrier: energy vector 'electricity' or 'heat'
        :param tech: 'pv','wt','owc','biomass_system','biomass_boiler'
        :param power: nominal power [kW]
        :param status: 'new' or 'old'
        :param user: list of ID user physically connected
        :param decay: production annual decay (for example 0.006)
        :param life_time: useful life of the plant [year]
        :param p_con: connections point with the network
        :param cost_kW: [€/kW]
        :param oem_cost_kW: maintenance costs [€/year/kW]
        :param inc: type of incentive for example bonus50%
        :param dur_inc:duration incentive [year]
        :param aux_components: list of auxiliary components
        :param mode:1 based on demand (BiomassBoiler is flexible); 2 based on biomass flow rate
        :param power_el_absorbed: list of [electrical power absorbed in working fase [kW],electrical power absorbed on power on [kW]] 
        :param eff: heating efficiency range:list of [Efficiency at reduced heating power,Efficiency at rated heating power] 
        :param fuel: wood pellets,wood chips,other
        :biomass_flow_rate_max: maximux admissible biomass flow rate [kg/h]
        :biomass_flow_rate_min: minimum admissible biomass flow rate [kg/h]
        :param pm_emissions_rate: particulate emissions [mg/Nm3]
        :param co_emissions_rate: CO emissions [mg/Nm3]
        
     
        """
        if self.category == 'flex':
            self.mode = 1
        else:
            self.mode = 2
        self.power_el_absorbed = power_el_absorbed
        self.eff = eff
        self.biomass_flow_rate_max = biomass_flow_rate_max
        self.biomass_flow_rate_min = biomass_flow_rate_min
        self.biomass = biomass
        self.pm_emissions_rate = pm_emmissiosn_rate
        self.co_emissions_rate = co_emissions_rate

        self.ignition_phase_duration = 0.16  # [h] about 10 minutes

    def efficiency_curve(self):
        """

        :return:efficiency vs power ratio
        """
        x = np.array(
            [self.biomass_flow_rate_min / conv_hour_sec * self.biomass.lhv * np.mean(np.array(self.eff)) / self.power,
             1])
        y = np.array(self.eff)
        f = interp1d(x, y, kind='linear')

        return f

    def controller(self, time, d=None, biomass_flow_rate=None):
        """

        :param d: heat demand [kW]
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :param biomass_flow_rate:[kg/time] (it is mandatory only in mode 2)
        :return:
        p: [kW]
        surplus: [kW]
        unmet: [kW]
        biomass_flow_rate: biomass consumption [kg/time]
        biomass_flow_rate_surplus: [kg/time]
        biomass_flow_rate_deficit: [kg/time]
        """

        power_max = self.biomass_flow_rate_max / conv_hour_sec * self.biomass.lhv * np.mean(np.array(self.eff))
        power_min = self.biomass_flow_rate_min / conv_hour_sec * self.biomass.lhv * np.mean(np.array(self.eff))
        if self.mode == 1:
            if d >= power_max:
                p = power_max
            elif d <= power_min and d != 0:
                p = power_min
            else:
                p = d

            if d==0:
                p=0
                biomass_flow_rate=0
            else:
                power_ratio = p / self.power
                f_eff = self.efficiency_curve()
                eff = f_eff(power_ratio)
                q_abs = p / eff
                biomass_flow_rate = q_abs * time * conv_hour_sec / self.biomass.lhv

            if d >= p:
                surplus = 0
                unmet = d - p
            else:
                surplus = p - d
                unmet = 0

            return p, surplus, unmet, biomass_flow_rate

        else:
            biomass_flow_rate_surplus = 0
            biomass_flow_rate_deficit = 0
            if biomass_flow_rate >= self.biomass_flow_rate_max * time:
                biomass_flow_rate_surplus = biomass_flow_rate - self.biomass_flow_rate_max * time
                biomass_flow_rate = self.biomass_flow_rate_max * time
                biomass_flow_rate_deficit = 0
            elif biomass_flow_rate <= self.biomass_flow_rate_min * time:
                biomass_flow_rate_deficit = self.biomass_flow_rate_min * time - biomass_flow_rate
                biomass_flow_rate = 0
                biomass_flow_rate_surplus = 0

            if biomass_flow_rate == 0:
                p = 0
            else:
                power_ratio = biomass_flow_rate / self.biomass_flow_rate_max
                f_eff = self.efficiency_curve()
                eff = f_eff(power_ratio)
                p = biomass_flow_rate / (time * conv_hour_sec) * self.biomass.lhv * eff

            return p, biomass_flow_rate, biomass_flow_rate_surplus, biomass_flow_rate_deficit

    def compute_output(self,  time,d=None, biomass_flow_rate=None):
        """

        :param d: heat demand [kW]
        :param time: 1 if hourly analysis, 0.25 if quarterly analysis
        :param biomass_flow_rate: [kg/time] (it is mandatory only in mode 2)
        :return:
        """

        if d is not None:
            len_ref=len(d)
        else:
            len_ref=len(biomass_flow_rate)


        p = np.zeros((len_ref,))
        surplus = np.zeros((len_ref,))
        unmet = np.zeros((len_ref,))
        biomass_cons = np.zeros((len_ref,))
        power_el_absorbed = np.zeros((len_ref,))
        tCO2 = np.zeros((len_ref,))
        co_emissions = np.zeros((len_ref,))
        pm_emissions = np.zeros((len_ref,))
        biomass_flow_rate_surplus = np.zeros((len_ref,))
        biomass_flow_rate_deficit = np.zeros((len_ref,))
        check = np.zeros((len_ref,))

        check_previous = 0
        for i in range(len_ref):

            if self.mode == 1:
                p_i, surplus_i, unmet_i, biomass_cons_i = self.controller(d=d[i], time=time)
                biomass_flow_rate_surplus_i = 0
                biomass_flow_rate_deficit_i = 0
            else:
                p_i, biomass_cons_i, biomass_flow_rate_surplus_i, biomass_flow_rate_deficit_i = self.controller( time=time,
                    biomass_flow_rate=biomass_flow_rate[i])
                surplus_i=None
                unmet_i=None

            tCO2_i = biomass_cons_i / 1000 * self.biomass.tCO2_emission

            if p_i == 0:
                check = 0
            else:
                check = 1

            if check_previous == 0 and check == 1:
                power_el_absorbed_i = self.power_el_absorbed[1]
            elif check == 0:
                power_el_absorbed_i = 0
            else:
                power_el_absorbed_i = self.power_el_absorbed[0]

            check_previous = check

            pm_emissions_i = 0
            co_emissions_i = 0

            p[i] = p_i
            surplus[i] = surplus_i
            unmet[i] = unmet_i
            biomass_cons[i] = biomass_cons_i
            power_el_absorbed[i] = power_el_absorbed_i
            tCO2[i] = tCO2_i
            co_emissions[i] = co_emissions_i
            pm_emissions[i] = pm_emissions_i
            biomass_flow_rate_surplus[i] = biomass_flow_rate_surplus_i
            biomass_flow_rate_deficit[i] = biomass_flow_rate_deficit_i

        self.en_perf_evolution['heat'] = {}
        self.en_perf_evolution['heat']['prod'] = p
        self.en_perf_evolution['heat']['surplus'] = surplus
        self.en_perf_evolution['heat']['unmet'] = unmet

        return p, surplus, unmet, biomass_cons, power_el_absorbed, tCO2, co_emissions, pm_emissions, biomass_flow_rate_surplus, biomass_flow_rate_deficit

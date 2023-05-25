from scipy import interpolate
import numpy as np


class Bess:
    def __init__(self, df_chargecurve, df_dischargecurve, t_ref, q_ref, v_charge_cut_off, v_cut_off, i_max_cell,
                 i_min_cell, n_series, n_parallel, soc_max, soc_min, eff, capacity, cost, oem_cost_kWh, replace,
                 replacement_year, plant):
        """

        :param df_chargecurve: charging curve on the manufacturer's spec sheet
        :param df_dischargecurve: discharge curve on the manufacturer's spec sheet
        :param t_ref: reference temperatures for which battery is tested [°C]
        :param q_ref: cell energy capacity [Ah]
        :param v_charge_cut_off:voltage at which battery is considered fully charged [V]
        :param v_cut_off: voltage at which battery is considered fully discharged [V]
        :param i_max_cell: max. current per cell charging and discharge [A]
        :param i_min_cell: min. current per cell charging and discharge [A]
        :param n_series: module in series
        :param n_parallel: module in parallel
        :param soc_max: max. state of charge
        :param soc_min: min. state of charge
        :param eff:charging/discharge efficiency
        :param capacity: [kWh]
        :param cost [€/kWh]
        :param oem_cost_kWh :[kWh/year]
        :param replace: (True or False)
        :param replacement_year
        :param plant: list of production systems physically connected
        """
        self.df_chargecurve = df_chargecurve
        self.df_dischargecurve = df_dischargecurve
        self.t_ref = t_ref
        self.q_ref = q_ref
        self.v_charge_cut_off = v_charge_cut_off
        self.v_cut_off = v_cut_off
        self.i_max_cell = i_max_cell
        self.i_min_cell = i_min_cell
        self.n_series = n_series
        self.n_parallel = n_parallel
        self.soc_max = soc_max
        self.soc_min = soc_min

        self.i_surplus_cell = 0
        self.q_ch = self.df_chargecurve['q']
        self.v_ch =self.df_chargecurve['v']
        self.q_dt0 = self.df_dischargecurve['q0']
        self.v_dt0 = self.df_dischargecurve['v0']
        self.q_dt1 = self.df_dischargecurve['q1']
        self.v_dt1 = self.df_dischargecurve['v1']
        self.q_dt2 = self.df_dischargecurve['q2']
        self.v_dt2 = self.df_dischargecurve['v2']
        self.q_dt3 = self.df_dischargecurve['q3']
        self.v_dt3 = self.df_dischargecurve['v3']
        self.q_dt4 = self.df_dischargecurve['q4']
        self.v_dt4 = self.df_dischargecurve['v4']
        self.f_charge = interpolate.interp1d(self.q_ch, self.v_ch, bounds_error=False, fill_value=self.v_charge_cut_off)
        self.f_t0 = interpolate.interp1d(self.q_dt0, self.v_dt0, bounds_error=False, fill_value=self.v_cut_off)
        self.f_t1 = interpolate.interp1d(self.q_dt1, self.v_dt1, bounds_error=False, fill_value=self.v_cut_off)
        self.f_t2 = interpolate.interp1d(self.q_dt2, self.v_dt2, bounds_error=False, fill_value=self.v_cut_off)
        self.f_t3 = interpolate.interp1d(self.q_dt3, self.v_dt3, bounds_error=False, fill_value=self.v_cut_off)
        self.f_t4 = interpolate.interp1d(self.q_dt4, self.v_dt4, bounds_error=False, fill_value=self.v_cut_off)
        self.eff = eff
        self.capacity = capacity
        self.cost = cost
        self.oem_cost_kWh = oem_cost_kWh
        self.replace = replace
        self.replacement_year = replacement_year
        self.plant = plant
        self.en_perf_evolution = {}  # {'soc':state of charge,'stored':stored power [kW],'supplied':supplied stored [kW]} level 0 and level 1

    def compute_cell(self, i_tot, t, soc_cell, time):
        """

        :param i_tot: current of charging(+)/discharge(-) [A]
        :param t:Temperature [°C]
        :param soc_cell: cell state of charge
        :param time: [s]
        :return:i_cell
                i_surplus_cell: excess current [A]
                i_deficit_cell: deficit current [A]
                soc_cell: new cell state of charge
                v_cell: cell voltage [V]
                q_discharge_cell_perc: reduction in state of charge [%]
                q_cell:cell state of charge [Ah]
                q_max_cell: maximum cell state of charge [Ah]
        """
        if t >= self.t_ref[0]:
            q_max_cell = self.q_ref[0]
        if self.t_ref[1] <= t < self.t_ref[0]:
            q_max_cell = ((t - self.t_ref[0]) * self.q_ref[1] + (self.t_ref[1] - t) * self.q_ref[0]) / (
                    self.t_ref[1] - self.t_ref[0])
        if self.t_ref[2] <= t < self.t_ref[1]:
            q_max_cell = ((t - self.t_ref[1]) * self.q_ref[2] + (self.t_ref[2] - t) * self.q_ref[1]) / (
                    self.t_ref[2] - self.t_ref[1])
        if self.t_ref[3] <= t < self.t_ref[2]:
            q_max_cell = ((t - self.t_ref[2]) * self.q_ref[3] + (self.t_ref[3] - t) * self.q_ref[2]) / (
                    self.t_ref[3] - self.t_ref[2])
        if t < self.t_ref[3]:
            q_max_cell = ((t - self.t_ref[3]) * self.q_ref[4] + (self.t_ref[4] - t) * self.q_ref[3]) / (
                    self.t_ref[4] - self.t_ref[3])

        i_cell = i_tot / self.n_parallel
        q_cell = soc_cell * q_max_cell
        v_cell = np.array([0, 0])
        i_deficit_cell = 0

        if i_cell > 0:  # charging
            if i_cell > self.i_max_cell:
                i_deficit_cell = 0
                i_surplus_cell = i_cell - self.i_max_cell
                i_cell = self.i_max_cell
            elif i_cell < self.i_min_cell:
                i_deficit_cell = 0
                i_surplus_cell = i_cell
                i_cell = 0

            else:
                i_deficit_cell = 0
                i_surplus_cell = 0

            q_max_cell = self.q_ref[0]
            soc_cell = soc_cell + i_cell * time / q_max_cell
            if soc_cell < self.soc_min:
                soc_deficit = soc_cell - self.soc_min
                soc_cell = self.soc_min
                i_defect = soc_deficit * q_max_cell
                i_deficit_cell = i_deficit_cell + i_defect
            if soc_cell > self.soc_max:
                soc_surplus = soc_cell - self.soc_max
                soc_cell = self.soc_max
                i_excess = soc_surplus * q_max_cell
                i_surplus_cell = i_surplus_cell + i_excess

            q_cell = soc_cell * q_max_cell
            v_cell = np.array([q_cell, self.f_charge(q_cell)])

        if i_cell <= 0:  # discharge or zero current
            if i_tot <= 0:
                if abs(i_cell) > self.i_max_cell:
                    i_surplus_cell = 0
                    i_deficit_cell = -(abs(i_cell) - self.i_max_cell)
                    i_cell = -self.i_max_cell
                elif abs(i_cell) < self.i_min_cell:
                    i_surplus_cell = 0
                    i_deficit_cell = i_cell
                    i_cell = 0
                else:
                    i_surplus_cell = 0
                    i_deficit_cell = 0
            else:
                i_surplus_cell = i_surplus_cell
                i_deficit_cell = 0

            soc_cell = soc_cell + i_cell * time / q_max_cell
            if soc_cell < self.soc_min:
                soc_deficit = soc_cell - self.soc_min
                soc_cell = self.soc_min
                i_defect = soc_deficit * q_max_cell
                i_deficit_cell = i_deficit_cell + i_defect
            if soc_cell > self.soc_max:
                soc_surplus = soc_cell - self.soc_max
                soc_cell = self.soc_max
                i_excess = soc_surplus * q_max_cell
                i_surplus_cell = i_surplus_cell + i_excess

            if t >= self.t_ref[0]:
                q_cell = (1 - soc_cell) * q_max_cell
                v_cell = np.array([q_cell, self.f_t0(q_cell)])
                q_cell = q_max_cell - q_cell
                if i_cell == 0:
                    v_cell = np.array([q_cell, ((self.f_charge(q_cell) + v_cell[1]) / 2)])

            if self.t_ref[1] <= t < self.t_ref[0]:
                q_cell = (1 - soc_cell) * q_max_cell
                v_sup = np.array([q_cell, self.f_t0(q_cell)])
                v_inf = np.array([q_cell, self.f_t1(q_cell)])
                v_cell = ((t - self.t_ref[0]) * v_inf + (self.t_ref[1] - t) * v_sup) / (self.t_ref[1] - self.t_ref[0])

                if q_cell > self.q_ref[1] and t != self.t_ref[1]:
                    volt_sup_i = np.array([q_cell, self.f_t0(self.q_ref[1])])
                    volt_inf_i = np.array([q_cell, self.f_t1(self.q_ref[1])])
                    volt_point_one_line = ((t - self.t_ref[0]) * volt_inf_i + (self.t_ref[1] - t) * volt_sup_i) / (
                            self.t_ref[1] - self.t_ref[0])
                    line_coeff = (self.v_cut_off - volt_point_one_line[1]) / (q_max_cell - self.q_ref[1])
                    v_cell = np.array([q_cell, line_coeff * (q_cell - self.q_ref[1]) + volt_point_one_line[1]])
                q_cell = q_max_cell - q_cell

                if i_cell == 0:
                    v_cell = np.array([q_cell, (self.f_charge(q_cell) + v_cell[1]) / 2])

            if self.t_ref[2] <= t < self.t_ref[1]:
                q_cell = (1 - soc_cell) * q_max_cell
                v_sup = np.array([q_cell, self.f_t1(q_cell)])
                v_inf = np.array([q_cell, self.f_t2(q_cell)])
                v_cell = ((t - self.t_ref[1]) * v_inf + (self.t_ref[2] - t) * v_sup) / (self.t_ref[2] - self.t_ref[1])

                if q_cell > self.q_ref[2] and t != self.t_ref[2]:
                    volt_sup_i = np.array([q_cell, self.f_t1(self.q_ref[2])])
                    volt_inf_i = np.array([q_cell, self.f_t2(self.q_ref[2])])
                    volt_point_one_line = ((t - self.t_ref[1]) * volt_inf_i + (self.t_ref[2] - t) * volt_sup_i) / (
                            self.t_ref[2] - self.t_ref[1])
                    line_coeff = (self.v_cut_off - volt_point_one_line[1]) / (q_max_cell - self.q_ref[2])
                    v_cell = np.array([q_cell, line_coeff * (q_cell - self.q_ref[2]) + volt_point_one_line[1]])
                q_cell = q_max_cell - q_cell

                if i_cell == 0:
                    v_cell = np.array([q_cell, (self.f_charge(q_cell) + v_cell[1]) / 2])

            if self.t_ref[3] <= t < self.t_ref[2]:
                q_cell = (1 - soc_cell) * q_max_cell
                v_sup = np.array([q_cell, self.f_t2(q_cell)])
                v_inf = np.array([q_cell, self.f_t3(q_cell)])
                v_cell = ((t - self.t_ref[2]) * v_inf + (self.t_ref[3] - t) * v_sup) / (self.t_ref[3] - self.t_ref[2])

                if q_cell > self.q_ref[3] and t != self.t_ref[3]:
                    volt_sup_i = np.array([q_cell, self.f_t2(self.q_ref[3])])
                    volt_inf_i = np.array([q_cell, self.f_t3(self.q_ref[3])])
                    volt_point_one_line = ((t - self.t_ref[2]) * volt_inf_i + (self.t_ref[3] - t) * volt_sup_i) / (
                            self.t_ref[3] - self.t_ref[2])
                    line_coeff = (self.v_cut_off - volt_point_one_line[1]) / (q_max_cell - self.q_ref[3])
                    v_cell = np.array([q_cell, line_coeff * (q_cell - self.q_ref[3]) + volt_point_one_line[1]])
                q_cell = q_max_cell - q_cell

                if i_cell == 0:
                    v_cell = np.array([q_cell, (self.f_charge(q_cell) + v_cell[1]) / 2])

            if t < self.t_ref[3]:
                q_cell = (1 - soc_cell) * q_max_cell
                v_sup = np.array([q_cell, self.f_t3(q_cell)])
                v_inf = np.array([q_cell, self.f_t4(q_cell)])
                v_cell = ((t - self.t_ref[3]) * v_inf + (self.t_ref[4] - t) * v_sup) / (self.t_ref[4] - self.t_ref[3])

                if q_cell > self.q_ref[4] and t != self.t_ref[4]:
                    v_sup = np.array([q_cell, self.f_t4(self.q_ref[4])])
                    v_inf = np.array([q_cell, self.f_t3(self.q_ref[4])])
                    volt_point_one_line = ((t - self.t_ref[3]) * v_inf + (self.t_ref[4] - t) * v_sup) / (
                            self.t_ref[4] - self.t_ref[3])
                    line_coeff = (self.v_cut_off - volt_point_one_line[1]) / (q_max_cell - self.q_ref[4])
                    v_cell = np.array([q_cell, line_coeff * (q_cell - self.q_ref[4]) + volt_point_one_line[1]])
                q_cell = q_max_cell - q_cell
                if i_cell == 0:
                    volt_discharge_i = v_cell[1]
                    v_cell = np.array([q_cell, (self.f_charge(q_cell) + volt_discharge_i) / 2])

        q_discharge_cell_perc = (q_max_cell - q_cell) * 100 / (self.q_ref[0])
        return i_cell, i_surplus_cell, i_deficit_cell, soc_cell, v_cell[1], q_discharge_cell_perc, q_cell, q_max_cell

    def compute_module(self, i_tot, t, soc_cell, time):
        '''

        :param i_tot: current of charging(+)/discharge(-) [A]
        :param t:Temperature [°C]
        :param soc_cell: cell state of charge
        :param time: [s]
        :return:v_tot: battery voltage [V]
                q_tot: battery state of charge [Ah]
                q_discharge_tot: reduction in state of charge [Ah]
                i_surplus: battery excess current [A]
                i_deficit: battery deficit current [A]
        '''

        i_cell, i_surplus_cell, i_deficit_cell, soc_cell, v_cell, q_discharge_cell_perc, q_cell, q_max_cell = self.compute_cell(
            i_tot=i_tot, t=t, soc_cell=soc_cell, time=time)
        v_tot = v_cell * self.n_series
        q_tot = q_cell * self.n_parallel
        q_discharge_tot = (q_max_cell - q_cell) * self.n_parallel
        i_surplus = i_surplus_cell * self.n_parallel
        i_deficit = i_deficit_cell * self.n_parallel

        return v_tot, q_tot, q_discharge_tot, i_surplus, i_deficit

    def controller(self, surplus, deficit, v_tot, soc, q_final_tot, time):
        """

        :param surplus: input power to battery [W]
        :param deficit: power required from battery [W]
        :param v_tot: battery voltage [V]
        :param soc: battery state of charge
        :param q_final_tot: battery state of charge [Ah]
        :param time: [s]
        :return:i_storage: current of charging(+)/discharge(-) [A]
                e_supply: supllied energy [Wh]
                e_stored: stored energy [Wh]
                e_deficit: deficit energy from battery[Wh]
                e_surplus: excess energy from battery [Wh]
                e_lost: energy lost during charging/discharge [Wh]
        """
        i_max_tot = self.i_max_cell * self.n_parallel
        i_min_tot = self.i_min_cell * self.n_parallel
        surplus = surplus * time
        deficit = deficit * time
        if surplus > 0:
            surplus_b = surplus * self.eff
            if soc < self.soc_max:
                a_storage = (self.soc_max - soc) * q_final_tot * v_tot  # [Wh]
                if a_storage <= surplus_b:
                    e_supply = 0
                    e_deficit = 0
                    i_availability = a_storage / (v_tot * time)
                    i_storage = min(i_availability, i_max_tot)
                    if i_storage < i_min_tot:
                        i_storage = 0
                    e_stored = i_storage * v_tot * time / self.eff
                    e_lost = e_stored * (1 - self.eff)
                    e_surplus = surplus - e_stored

                else:
                    e_supply = 0
                    e_deficit = 0
                    i_availability = surplus * self.eff / (v_tot * time)
                    i_storage = min(i_availability, i_max_tot)
                    if i_storage < i_min_tot:
                        i_storage = 0
                    e_stored = i_storage * v_tot * time / self.eff
                    e_lost = e_stored * (1 - self.eff)
                    e_surplus = surplus - e_stored

            else:
                e_stored = 0
                e_lost = 0
                e_surplus = surplus
                e_supply = 0
                e_deficit = 0
                i_storage = 0

        elif deficit > 0:
            if soc > self.soc_min:
                supply = (soc - self.soc_min) * q_final_tot * v_tot  # [Wh]
                if supply > deficit / self.eff:
                    e_supply = deficit  # potenzialmente fornibile
                    e_stored = 0
                    e_surplus = 0
                    i_supply = e_supply / (self.eff * v_tot * time)
                    i_storage = min(i_supply, i_max_tot)
                    if i_storage < i_min_tot:
                        i_storage = 0
                    e_supply = i_storage * v_tot * self.eff * time
                    e_lost = i_storage * v_tot * (1 - self.eff) * time
                    e_deficit = deficit - e_supply
                    i_storage = -i_storage
                else:
                    e_supply = supply
                    e_stored = 0
                    e_surplus = 0
                    i_supply = e_supply / (v_tot * time)
                    i_storage = min(i_supply, i_max_tot)
                    if i_storage < i_min_tot:
                        i_storage = 0
                    e_supply = i_storage * v_tot * self.eff * time
                    e_lost = i_storage * v_tot * (1 - self.eff) * time
                    e_deficit = deficit - e_supply
                    i_storage = - i_storage
            else:
                e_supply = 0
                e_lost = 0
                e_deficit = deficit
                e_stored = 0
                e_surplus = 0
                i_storage = 0
        else:
            e_supply = 0
            e_lost = 0
            e_deficit = 0
            e_stored = 0
            e_surplus = 0
            i_storage = 0

        return i_storage, e_supply, e_stored, e_deficit, e_surplus, e_lost

    def compute_output(self, surplus, deficit, t, soc_cell_in, time, i_tot_in=0, t_in=25):

        i_cell_i, i_surplus_cell_i, i_deficit_cell_i, soc_cell_i, v_cell_i, q_discharge_cell_perc_i, q_cell_i, q_max_cell_i = self.compute_cell(
            i_tot=i_tot_in, t=t_in, soc_cell=soc_cell_in, time=time)
        v_tot_i, q_tot_i, q_discharge_tot_i, i_surplus_i, i_deficit_i = self.compute_module(i_tot=i_tot_in, t=t_in,
                                                                                            soc_cell=soc_cell_in,
                                                                                            time=time)
        i_storage = np.zeros((len(surplus),))
        e_supply = np.zeros((len(surplus),))
        e_stored = np.zeros((len(surplus),))
        e_deficit = np.zeros((len(surplus),))
        e_surplus = np.zeros((len(surplus),))
        e_lost = np.zeros((len(surplus),))
        v_tot = np.zeros((len(surplus),))
        q_tot = np.zeros((len(surplus),))
        q_discharge_tot = np.zeros((len(surplus),))
        i_surplus = np.zeros((len(surplus),))
        i_deficit = np.zeros((len(surplus),))

        for i in range(len(surplus)):
            i_storage_i, e_supply_i, e_stored_i, e_deficit_i, e_surplus_i, e_lost_i = self.controller(
                surplus=surplus[i] * 1000, deficit=deficit[i] * 1000, v_tot=v_tot_i, soc=soc_cell_i,
                q_final_tot=q_max_cell_i * self.n_parallel, time=time)
            i_cell_i, i_surplus_cell_i, i_deficit_cell_i, soc_cell_i, v_cell_i, q_discharge_cell_perc_i, q_cell_i, q_max_cell_i = self.compute_cell(
                i_tot=i_storage_i, t=t[i], soc_cell=soc_cell_i, time=time)
            v_tot_i, q_tot_i, q_discharge_tot_i, i_surplus_i, i_deficit_i = self.compute_module(i_tot=i_storage_i,
                                                                                                t=t[i],
                                                                                                soc_cell=soc_cell_i,
                                                                                                time=time)
            i_storage[i], e_supply[i], e_stored[i], e_deficit[i], e_surplus[i], e_lost[
                i] = i_storage_i, e_supply_i, e_stored_i, e_deficit_i, e_surplus_i, e_lost_i
            v_tot[i], q_tot[i], q_discharge_tot[i], i_surplus[i], i_deficit[
                i] = v_tot_i, q_tot_i, q_discharge_tot_i, i_surplus_i, i_deficit_i

        return i_storage, e_supply, e_stored, e_deficit, e_surplus, e_lost, v_tot, q_tot, q_discharge_tot, i_surplus, i_deficit

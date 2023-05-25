import numpy as np
import numpy_financial as npf


def order_by_kW(plant):
    parameter = plant.power
    return parameter


class Prosumer:
    def __init__(self, id, plant, user, list_carrier):
        """
        to simulate user physically connected to power plants.
        :param id: identification code.
        :param plant: list of distributed generation plants.
        :param user: list of user physically connected to plants.
        :param list_carrier: list of carrier ['electricity','heat'] #il vettore principale deve essere messo prima
        """
        self.id = id
        self.user = user
        self.plant = plant
        self.list_carrier = list_carrier

        self.en_perf_evolution = {}  # {p_tot, d_tot, self_cons, surplus, deficit} level 0 plus {surplus_1, deficit_1} level 1
        self.ec_perf_evolution = {}  # [NPV,cash_flow,cum_cash_flow,pbp]

    def compute_list_aux_components(self):
        list_aux_components = []
        for plant in self.plant:
            if plant.aux_components:
                list_aux_components += plant.aux_components
            else:
                list_aux_components = list_aux_components
        self.list_aux_components = list_aux_components
        return list_aux_components

    def energy_perfomance(self, time):
        """

        :return:production [kW]
        :return:demand [kW]
        :return:self_cons [kW]
        :return: surplus [kW]
        :return: deficit [kW]

        """
        if self.list_carrier[0] == 'electricity':
            p_heat_cog = 0
            for carrier in self.list_carrier:
                d_tot = 0
                for user in self.user:
                    if carrier in user.carrier:
                        d_tot += user.en_perf_evolution[carrier]
                    else:
                        d_tot += 0

                # if carrier == 'heat':
                #     if type(p_heat_cog) != int:
                #         p_inflex = p_heat_cog
                #     else:
                #         p_inflex = np.zeros((len(d_tot),))
                # else:
                #     p_inflex = np.zeros((len(d_tot),))

                p_inflex=np.zeros((len(d_tot),))
                p_flex = np.zeros((len(d_tot),))
                p_tot = np.zeros((len(d_tot),))

                plant_inflex = []
                plant_flex = []
                for plant in self.plant:
                    if carrier in plant.carrier:
                        if plant.category == 'inflex':
                            plant_inflex.append(plant)
                        else:
                            plant_flex.append(plant)
                    if carrier == 'heat':
                        if plant.carrier == 'electricity' and plant.cogeneration == True:
                            plant_inflex.append(plant)

                if plant_inflex:
                    for plant in plant_inflex:
                        p_inflex += plant.en_perf_evolution[carrier]['prod']
                        if carrier == 'electricity' and plant.cogeneration is True:
                            p_heat_cog += plant.en_perf_evolution['heat']['prod']

                self_cons = np.zeros((len(d_tot),))
                surplus = np.zeros((len(d_tot),))
                unmet = np.zeros((len(d_tot),))
                for dt in range(len(d_tot)):
                    p = p_inflex[dt]
                    d = d_tot[dt]
                    if p >= d:
                        surplus[dt] = p - d
                        unmet[dt] = 0

                    else:
                        surplus[dt] = 0
                        unmet[dt] = d - p

                if plant_flex:
                    plant_flex.sort(key=order_by_kW)
                    for item in range(len(plant_flex)):
                        output = plant_flex[item].compute_output(time=time, d=unmet)
                        p_flex += plant_flex[item].en_perf_evolution[carrier]['prod']
                        surplus += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        unmet = plant_flex[item].en_perf_evolution[carrier]['unmet']
                        if carrier == 'electricity' and plant_flex[item].cogeneration is True:
                            p_heat_cog += plant_flex[item].en_perf_evolution['heat']['prod']

                if type(p_heat_cog) != int:
                    self.en_perf_evolution['heat_cog'] = {}
                    self.en_perf_evolution['heat_cog']['prod'] = p_heat_cog
                    if carrier == 'heat':
                        self_cons_heat = np.zeros((len(d_tot),))
                        for dt in range(len(d_tot)):
                            p = p_heat_cog[dt]
                            d = d_tot[dt]
                            if p >= d:
                                self_cons_heat[dt] = d
                            else:
                                self_cons_heat[dt] = p
                        self.en_perf_evolution['heat_cog']['self_cons'] = self_cons_heat

                for dt in range(len(d_tot)):
                    p = p_inflex[dt] + p_flex[dt]
                    d = d_tot[dt]
                    p_tot[dt] = p
                    if p >= d:
                        self_cons[dt] = d

                    else:
                        self_cons[dt] = p

                self.en_perf_evolution[carrier] = {}
                self.en_perf_evolution[carrier]['prod'] = p_tot
                self.en_perf_evolution[carrier]['dem'] = d_tot
                self.en_perf_evolution[carrier]['self_cons'] = self_cons
                self.en_perf_evolution[carrier]['surplus'] = surplus
                self.en_perf_evolution[carrier]['unmet'] = unmet
                self.en_perf_evolution[carrier]['prod_inflex'] = p_inflex
                self.en_perf_evolution[carrier]['prod_flex'] = p_flex
        else:
            p_el_cog = 0
            for carrier in self.list_carrier:
                d_tot = 0
                for user in self.user:
                    if carrier in user.carrier:
                        d_tot += user.en_perf_evolution[carrier]
                    else:
                        d_tot += 0
                #
                # if carrier == 'electricity':
                #     if type(p_el_cog) != int:
                #         p_inflex = p_el_cog
                #     else:
                #         p_inflex = np.zeros((len(d_tot),))
                # else:
                #     p_inflex = np.zeros((len(d_tot),))

                p_inflex = np.zeros((len(d_tot),))

                p_flex = np.zeros((len(d_tot),))
                p_tot = np.zeros((len(d_tot),))

                plant_inflex = []
                plant_flex = []
                for plant in self.plant:
                    if carrier in plant.carrier:
                        if plant.category == 'inflex':
                            plant_inflex.append(plant)
                        else:
                            plant_flex.append(plant)
                    if carrier == 'electricity':
                        if plant.carrier == 'heat' and plant.cogeneration == True:
                            plant_inflex.append(plant)
                if plant_inflex:
                    for plant in plant_inflex:
                        p_inflex += plant.en_perf_evolution[carrier]['prod']
                        if carrier == 'heat' and plant.cogeneration is True:
                            p_el_cog += plant.en_perf_evolution['electricity']['prod']

                self_cons = np.zeros((len(d_tot),))
                surplus = np.zeros((len(d_tot),))
                unmet = np.zeros((len(d_tot),))
                for dt in range(len(d_tot)):
                    p = p_inflex[dt]
                    d = d_tot[dt]
                    if p >= d:
                        surplus[dt] = p - d
                        unmet[dt] = 0

                    else:
                        surplus[dt] = 0
                        unmet[dt] = d - p

                if plant_flex:
                    plant_flex.sort(key=order_by_kW)
                    for item in range(len(plant_flex)):
                        output = plant_flex[item].compute_output(time=time, d=unmet)
                        p_flex += plant_flex[item].en_perf_evolution[carrier]['prod']
                        surplus += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        unmet = plant_flex[item].en_perf_evolution[carrier]['unmet']
                        if carrier == 'electricity' and plant_flex[item].cogeneration is True:
                            p_el_cog += plant_flex[item].en_perf_evolution['heat']['prod']

                if type(p_el_cog) != int:
                    self.en_perf_evolution['el_cog'] = {}
                    self.en_perf_evolution['el_cog']['prod'] = p_el_cog

                for dt in range(len(d_tot)):
                    p = p_inflex[dt] + p_flex[dt]
                    d = d_tot[dt]
                    p_tot[dt] = p
                    if p >= d:
                        self_cons[dt] = d

                    else:
                        self_cons[dt] = p

                self.en_perf_evolution[carrier] = {}
                self.en_perf_evolution[carrier]['prod'] = p_tot
                self.en_perf_evolution[carrier]['dem'] = d_tot
                self.en_perf_evolution[carrier]['self_cons'] = self_cons
                self.en_perf_evolution[carrier]['surplus'] = surplus
                self.en_perf_evolution[carrier]['unmet'] = unmet
                self.en_perf_evolution[carrier]['prod_inflex'] = p_inflex
                self.en_perf_evolution[carrier]['prod_flex'] = p_flex

        return self.en_perf_evolution

    def economic_perfomance(self, time, t_inv, down_payment_percentual, t_res, int_rate, pr_import, pr_export, tax,
                            value_cb=None):
        """

        :param time:1 if hourly analysis [h], 0.25 if quarterly analysis [1/4 h]
        :param t_inv: investment time horizon [year]
        :param down_payment_percentual: percentual of investment paid to zero year [%]
        :param t_res:residual time to pay off debt [year]
        :param int_rate:interest rate
        :param pr_import: purchase price of imported energy per vector [€/MWh] (dict)
        :param pr_export: sale price of exported energy per vector [€/MWh] (dict)
        :param tax: percentual of taxes on sold energy  [%]
        :return:
            cashflow
            cashflow_cum
            pbp: pay back period [year]
            r1:revenue from sold electricity [€/year]
            r2: revenue from savings on energy bill [€/year]
            r3: revenue from incentive on kWh [€/year]
            r4: revenue from incentive [€/year]
            c1: cost of resources [€/year]
            c2: cost of oem [€/year] (includes auxiliary components replacement)
            c3: tax to network system [€/year]
            c4:tax on energy sold [€/year]
            c5:depreciation (if all costs are not paid to year zero) [€/year]
        """

        energy_perfomance = self.energy_perfomance(time=time)

        list_carrier = self.list_carrier

        list_aux_components = self.compute_list_aux_components()
        plant_heat = []
        plant_el = []
        decay = {}
        for plant in self.plant:
            if 'electricity' is plant.carrier:
                plant_el.append(plant)
                if plant.cogeneration == True:
                    plant_heat.append(plant)
            elif 'heat' is plant.carrier:
                plant_heat.append(plant)
                if plant.cogeneration == True:
                    plant_el.append(plant)
        if plant_el:
            decay_aux_el = []
            for plant in plant_el:
                decay_aux_el.append(plant.decay)
            decay['electricity'] = np.mean(decay_aux_el)
        if plant_heat:
            decay_aux_heat = []
            for plant in plant_heat:
                decay_aux_heat.append(plant.decay)
            decay['heat'] = np.mean(decay_aux_heat)

        outflow = np.zeros((t_inv + 1,))
        inflow = np.zeros((t_inv + 1,))
        cashflow = np.zeros((t_inv + 1,))
        cashflow_cum = np.zeros((t_inv + 1,))
        r1 = np.zeros((t_inv + 1,))
        r2 = np.zeros((t_inv + 1,))
        r3 = np.zeros((t_inv + 1,))
        r4 = np.zeros((t_inv + 1,))
        c1 = np.zeros((t_inv + 1,))
        c2 = np.zeros((t_inv + 1,))
        c3 = np.zeros((t_inv + 1,))
        c4 = np.zeros((t_inv + 1,))
        c5 = np.zeros((t_inv + 1,))

        outflow0 = 0
        investment_cost = 0
        oem_cost = 0
        inc_kW = 0
        inc = 0
        gse_tax = 0
        for plant in self.plant:
            investment_cost += plant.cost
            inc_kW += plant.inc_kW
            outflow0 += (plant.cost - inc_kW) * down_payment_percentual / 100
            oem_cost += plant.oem_cost
            gse_tax += plant.gse_tax

        outflow[0] = outflow0
        inflow[0] = 0
        cashflow[0] = inflow[0] - outflow[0]
        cashflow_cum[0] = cashflow[0]
        r1[0] = 0
        r2[0] = 0
        r3[0] = 0
        c1[0] = 0
        c2[0] = 0
        c3[0] = 0
        c4[0] = 0
        c5[0] = 0

        for year in range(1, t_inv + 1):
            r1_i = 0
            r2_i = 0
            r3_i = 0
            r4_i = 0
            c1_i = 0
            c2_i = oem_cost
            for carrier in list_carrier:
                r1_i += sum(energy_perfomance[carrier]['surplus']) * time / (1000) * pr_export[carrier] * (
                        1 - decay[carrier]) ** (year - 1)
                r2_i += sum(energy_perfomance[carrier]['self_cons']) * time / (1000) * pr_import[carrier] * (
                        1 - decay[carrier]) ** (year - 1)

            for plant in self.plant:
                if year <= plant.dur_inc_kWh:
                    r3_i += sum(plant.en_perf_evolution[plant.carrier]['prod']) * time * plant.inc_kWh * (
                            1 - plant.decay) ** (
                                    year - 1)
                else:
                    r3_i += 0

                if year <= plant.dur_inc:
                    r4_i += plant.inc

                else:
                    r4_i += 0

                c1_i += plant.cost_resource
                c2_i += sum(plant.en_perf_evolution[plant.carrier]['prod']) * time * plant.oem_cost_kWh * (
                        1 - plant.decay) ** (
                                year - 1)

                if plant.cogeneration == True:
                    if plant.carrier == 'heat':
                        r4_i += plant.check_eff_calculate_cb(time=time, value_cb=value_cb,
                                                             self_cons_heat=None)[3] * (1 - plant.decay) ** (
                                        year - 1)
                    else:
                        r4_i += plant.check_eff_calculate_cb(time=time, value_cb=value_cb,
                                                             self_cons_heat=self.en_perf_evolution['heat_cog'][
                                                                 'self_cons'])[3] * (1 - plant.decay) ** (
                                        year - 1)

            for aux_components in list_aux_components:
                if aux_components.replacement_year == year:
                    c2_i += aux_components.replacement_cost
                else:
                    c2_i += 0

            c3_i = gse_tax
            c4_i = r1_i * tax / 100

            if down_payment_percentual != 100:
                if year <= t_res:
                    c5_i = investment_cost * (1 - down_payment_percentual / 100) / t_res
                else:
                    c5_i = 0
            else:
                c5_i = 0

            r1[year] = r1_i
            r2[year] = r2_i
            r3[year] = r3_i
            r4[year] = r4_i
            c1[year] = c1_i
            c2[year] = c2_i
            c3[year] = c3_i
            c4[year] = c4_i
            c5[year] = c5_i

            outflow[year] = c1_i + c2_i + c3_i + c4_i
            inflow[year] = r1_i + r2_i + r3_i + r4_i
            cashflow[year] = inflow[year] - outflow[year]
            cashflow_cum[year] = cashflow_cum[year - 1] + cashflow[year]

        NPV = npf.npv(int_rate, cashflow)
        pbp = outflow0 / np.mean(cashflow[1:])

        self.ec_perf_evolution['NPV'] = NPV
        self.ec_perf_evolution['cashflow'] = cashflow
        self.ec_perf_evolution['cashflow_cum'] = cashflow_cum
        self.ec_perf_evolution['pbp'] = pbp
        self.ec_perf_evolution['r1'] = r1
        self.ec_perf_evolution['r2'] = r2
        self.ec_perf_evolution['r3'] = r3
        self.ec_perf_evolution['r4'] = r4
        self.ec_perf_evolution['c1'] = c1
        self.ec_perf_evolution['c2'] = c2
        self.ec_perf_evolution['c3'] = c3
        self.ec_perf_evolution['c4'] = c4
        self.ec_perf_evolution['c5'] = c5

        return NPV, cashflow, cashflow_cum, pbp, r1, r2, r3, r4, c1, c2, c3, c4, c5

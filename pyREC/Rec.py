import numpy as np
import numpy_financial as npf


def order_by_kW(plant):
    parameter = plant.power
    return parameter


class Rec:
    def __init__(self, id, prosumer, consumer, rec_plant,list_carrier):
        """

        To simulate REC with consumers and prosumers
        :param id: identification code
        :param prosumer: list of prosumer
        :param consumer: list of consumer
        :param rec_plant: list of shared plants in the community
        :param list_carrier: lista dei vettori della comunità
        """

        self.id = id
        self.prosumer = prosumer
        self.consumer = consumer
        self.rec_plant = rec_plant
        self.list_carrier = list_carrier
        self.en_perf_evolution = {}
        self.ec_perf_evolution = {}

    def compute_list_aux_components(self):
        list_aux_components = []
        if self.rec_plant:
            for plant in self.rec_plant:
                if plant.aux_components:
                    list_aux_components += plant.aux_components
                else:
                    list_aux_components = list_aux_components
        self.list_aux_components = list_aux_components
        return list_aux_components

    def compute_members(self):
        if self.prosumer:
            n_prosumer = len(self.prosumer)
        else:
            n_prosumer = 0
        if self.consumer:
            n_consumer = len(self.consumer)
        else:
            n_consumer = 0

        n_member = n_prosumer + n_consumer

        return n_member, n_prosumer, n_consumer

    def energy_perfomance(self, time):

        # for prosumer in self.prosumer:
        #     energy_performance = prosumer.energy_perfomance(time=time)
        if self.list_carrier[0]=='electricity':
            p_heat_cog = 0
            for carrier in self.list_carrier:
                p_surplus = 0
                d_unmet = 0
                d_consumer = 0

                for prosumer in self.prosumer:
                    if carrier in prosumer.list_carrier:
                        p_surplus += prosumer.en_perf_evolution[carrier]['surplus']
                        d_unmet += prosumer.en_perf_evolution[carrier]['unmet']
                for consumer in self.consumer:
                    if carrier in consumer.carrier:
                        d_consumer += consumer.en_perf_evolution[carrier]

                if time==1:
                    len_ref = 8760
                else:
                    len_ref=35040


                if type(p_surplus)==int:
                    p_surplus = np.zeros((len_ref,))
                if type(d_consumer) == int:
                    d_consumer = np.zeros((len_ref,))
                if type(d_unmet) == int:
                    d_unmet = np.zeros((len_ref,))



                p_flex=np.zeros((len_ref,))
                p_inflex = np.zeros((len_ref,))

                plant_inflex = []
                plant_flex = []
                if self.rec_plant:
                    for plant in self.rec_plant:
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

                prod_tot = p_surplus + p_inflex
                prod_rec=p_inflex
                dem_tot = d_unmet + d_consumer
                surplus_tot = np.zeros((len_ref,))
                surplus_rec = np.zeros((len_ref,))
                unmet_tot = np.zeros((len_ref,))
                shared = np.zeros((len_ref,))
                for dt in range(len_ref):
                    p = prod_tot[dt]
                    p_rec=prod_rec[dt]
                    d = dem_tot[dt]
                    if p >= d:
                        surplus_tot[dt] = p - d
                        unmet_tot[dt] = 0
                    else:
                        surplus_tot[dt] = 0
                        unmet_tot[dt] = d - p
                    if p_rec>=d:
                        surplus_rec[dt]=p_rec-d
                    else:
                        surplus_rec[dt]=0

                if plant_flex:
                    plant_flex.sort(key=order_by_kW)
                    for item in range(len(plant_flex)):
                        output = plant_flex[item].compute_output(time=time, d=unmet_tot)
                        p_flex += plant_flex[item].en_perf_evolution[carrier]['prod']
                        surplus_tot += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        surplus_rec += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        unmet_tot = plant_flex[item].en_perf_evolution[carrier]['unmet']
                        if carrier == 'electricity' and plant_flex[item].cogeneration is True:
                            p_heat_cog += plant_flex[item].en_perf_evolution['heat']['prod']

                    prod_tot += p_flex
                    prod_rec+=p_rec

                for dt in range(len_ref):
                    p = prod_tot[dt]
                    d = dem_tot[dt]
                    if p >= d:
                        shared[dt] = d

                    else:
                        shared[dt] = p

                self.en_perf_evolution[carrier] = {}
                self.en_perf_evolution[carrier]['prod'] = prod_tot
                self.en_perf_evolution[carrier]['dem'] = dem_tot
                self.en_perf_evolution[carrier]['shared'] = shared
                self.en_perf_evolution[carrier]['surplus'] = surplus_tot
                self.en_perf_evolution[carrier]['surplus_rec'] = surplus_rec
                self.en_perf_evolution[carrier]['unmet'] = unmet_tot
                self.en_perf_evolution[carrier]['prod_inflex'] = p_inflex
                self.en_perf_evolution[carrier]['prod_flex'] = p_flex

                if type(p_heat_cog) != int:
                    self.en_perf_evolution['heat_cog'] = {}
                    self.en_perf_evolution['heat_cog']['prod'] = p_heat_cog
        else:
            p_el_cog = 0
            for carrier in self.list_carrier:
                p_surplus = 0
                d_unmet = 0
                d_consumer = 0
                p_flex = 0
                for prosumer in self.prosumer:
                    if carrier in prosumer.list_carrier:
                        p_surplus += prosumer.en_perf_evolution[carrier]['surplus']
                        d_unmet += prosumer.en_perf_evolution[carrier]['unmet']
                for consumer in self.consumer:
                    if carrier in consumer.carrier:
                        d_consumer += consumer.en_perf_evolution[carrier]

                len_ref = len(d_consumer)

                p_inflex = np.zeros((len_ref,))

                plant_inflex = []
                plant_flex = []
                if self.rec_plant:
                    for plant in self.rec_plant:
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

                prod_tot = p_surplus + p_inflex
                prod_rec = p_inflex
                dem_tot = d_unmet + d_consumer
                surplus_tot = np.zeros((len_ref,))
                surplus_rec = np.zeros((len_ref,))
                unmet_tot = np.zeros((len_ref,))
                shared = np.zeros((len_ref,))
                for dt in range(len_ref):
                    p = prod_tot[dt]
                    p_rec = prod_rec[dt]
                    d = dem_tot[dt]
                    if p >= d:
                        surplus_tot[dt] = p - d
                        unmet_tot[dt] = 0
                    else:
                        surplus_tot[dt] = 0
                        unmet_tot[dt] = d - p
                    if p_rec >= d:
                        surplus_rec[dt] = p_rec - d
                    else:
                        surplus_rec[dt] = 0

                if plant_flex:
                    plant_flex.sort(key=order_by_kW)
                    for item in range(len(plant_flex)):
                        output = plant_flex[item].compute_output(time=time, d=unmet_tot)
                        p_flex += plant_flex[item].en_perf_evolution[carrier]['prod']
                        surplus_tot += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        surplus_rec += plant_flex[item].en_perf_evolution[carrier]['surplus']
                        unmet_tot = plant_flex[item].en_perf_evolution[carrier]['unmet']
                        if carrier == 'heat' and plant_flex[item].cogeneration is True:
                            p_el_cog += plant_flex[item].en_perf_evolution['electricity']['prod']

                    prod_tot += p_flex
                    prod_rec += p_rec

                for dt in range(len_ref):
                    p = prod_tot[dt]
                    d = dem_tot[dt]
                    if p >= d:
                        shared[dt] = d

                    else:
                        shared[dt] = p

                self.en_perf_evolution[carrier] = {}
                self.en_perf_evolution[carrier]['prod'] = prod_tot
                self.en_perf_evolution[carrier]['dem'] = dem_tot
                self.en_perf_evolution[carrier]['shared'] = shared
                self.en_perf_evolution[carrier]['surplus'] = surplus_tot
                self.en_perf_evolution[carrier]['surplus_rec'] = surplus_rec
                self.en_perf_evolution[carrier]['unmet'] = unmet_tot
                self.en_perf_evolution[carrier]['prod_inflex'] = p_inflex
                self.en_perf_evolution[carrier]['prod_flex'] = p_flex

                if type(p_el_cog) != int:
                    self.en_perf_evolution['el_cog'] = {}
                    self.en_perf_evolution['el_cog']['prod'] = p_el_cog

        return self.en_perf_evolution

    def economic_perfomance(self, time, t_inv, down_payment_percentual, t_res, int_rate, inc_shared, pr_export, tax, p1,
                            p2):
        """

        :param time:1 if hourly analysis [h], 0.25 if quarterly analysis [1/4 h]
        :param t_inv: investment time horizon [year]
        :param down_payment_percentual: percentual of investment paid to zero year [%]
        :param t_res:residual time to pay off debt [year]
        :param int_rate:interest rate
        :param inc_shared: incentive of shared energy per vector [€/MWh] (dict)
        :param pr_export: sale price of exported energy per vector [€/MWh] (dict)
        :param tax: percentual of taxes on sold energy  [%]
        :param p1: partition coeff1 (until prosumer cashflow<0)
        :param p2: partition coeff2 (after prosumer cashflow>0)
        :return:
            cashflow
            cashflow_cum
            pbp: pay back period [year]
            r1:revenue from sold [€/year]
            r2: revenue from shared [€/year]
            r3: revenue from incentive on kWh produced [€/year]
            r4: revenue from incentive [€/year]
            c1: cost of resources [€/year]
            c2: cost of oem [€/year] (includes auxiliary components replacement)
            c3: tax to network system [€/year]
            c4: tax on energy sold [€/year]
            c5: depreciation (if all costs are not paid to year zero) [€/year]
            c6: cost of configuration REC by GSE [€/year]
            revenue_x_member: revenue for single member [€/member]
            revenue_x_member_x_year: revenue for single member per year [€/member/year]
        """

        energy_perfomance = self.energy_perfomance(time=time)


        list_aux_components = self.compute_list_aux_components()

        n_member, n_prosumer, n_consumer = self.compute_members()

        decay = {}
        for carrier in self.list_carrier:
            decay_aux = []
            if self.rec_plant:
                for plant in self.rec_plant:
                    if carrier is plant.carrier:
                        decay_plant = plant.decay
                        decay_aux.append(decay_plant)
            for prosumer in self.prosumer:
                for plant in prosumer.plant:
                    if carrier is plant.carrier:
                        decay_plant = plant.decay
                        decay_aux.append(decay_plant)

            decay[carrier] = np.mean(decay_aux)

        plant_heat = []
        plant_el = []
        decay = {}
        for plant in self.rec_plant:
            if 'electricity' is plant.carrier:
                plant_el.append(plant)
                if plant.cogeneration == True:
                    plant_heat.append(plant)
            elif 'heat' is plant.carrier:
                plant_heat.append(plant)
                if plant.cogeneration == True:
                    plant_el.append(plant)
        for prosumer in self.prosumer:
            for plant in prosumer.plant:
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
        c6 = np.zeros((t_inv + 1,))

        outflow0 = 0
        investment_cost = 0
        oem_cost = 0
        inc_kW = 0
        inc = 0
        gse_tax = 0
        if self.rec_plant:
            for plant in self.rec_plant:
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
        c6[0] = 0

        for year in range(1, t_inv + 1):
            r1_i = 0
            r2_i = 0
            r3_i = 0
            r4_i = 0
            c1_i = 0
            for carrier in self.list_carrier:
                r1_i += sum(energy_perfomance[carrier]['surplus_rec']) * time / 1000 * pr_export[carrier] * (
                        1 - decay[carrier]) ** (year - 1)
                r2_i += sum(energy_perfomance[carrier]['shared']) * time / 1000 * inc_shared[carrier] * (
                        1 - decay[carrier]) ** (year - 1)

            if 'electricity' in self.list_carrier:
                c6_i = 4*n_member
            else:
                c6_i = 0

            if self.rec_plant:
                for plant in self.rec_plant:
                    if year <= plant.dur_inc_kWh:
                        r3_i += sum(plant.en_perf_evolution[plant.carrier]['prod']) * plant.inc_kWh * (
                                1 - plant.decay) ** (
                                        year - 1)
                    else:
                        r3_i += 0

                    if year <= plant.dur_inc:
                        r4_i += plant.inc
                    else:
                        r4_i += 0

                    c1_i += plant.cost_resource

                c2_i = oem_cost
                for aux_components in list_aux_components:
                    if aux_components.replacement_year == year:
                        c2_i += aux_components.replacement_cost
                    else:
                        c2_i += 0

                c3_i = gse_tax

                if down_payment_percentual != 100:
                    if year <= t_res:
                        c5_i = investment_cost * (1 - down_payment_percentual / 100) / t_res
                    else:
                        c5_i = 0
                else:
                    c5_i = 0

            else:
                r3_i = 0
                r4_i = 0
                c1_i = 0
                c2_i = 0
                c3_i = 0
                c5_i = 0

            c4_i = r1_i * tax / 100

            r1[year] = r1_i
            r2[year] = r2_i
            r3[year] = r3_i
            r4[year] = r4_i
            c1[year] = c1_i
            c2[year] = c2_i
            c3[year] = c3_i
            c4[year] = c4_i
            c5[year] = c5_i
            c6[year] = c6_i

            outflow[year] = c1_i + c2_i + c3_i + c4_i
            inflow[year] = r1_i + r2_i + r3_i + r4_i
            cashflow[year] = inflow[year] - outflow[year]
            cashflow_cum[year] = cashflow_cum[year - 1] + cashflow[year]

        if self.prosumer:
            for prosumer in self.prosumer:
                cashflow_prosumer = np.zeros((t_inv + 1,))
                cashflow_cum_prosumer = np.zeros((t_inv + 1,))
                cashflow_prosumer[0] = prosumer.ec_perf_evolution['cashflow'][0]
                cashflow_cum_prosumer[0] = prosumer.ec_perf_evolution['cashflow_cum'][0]
                for year in range(1, t_inv + 1):
                    if prosumer.ec_perf_evolution['cashflow_cum'][year] <= 0:
                        cashflow_prosumer[year] = prosumer.ec_perf_evolution['cashflow'][year] + (p1 * r2[year]) / (
                            n_prosumer)
                        cashflow[year] = cashflow[year] - (p1 * r2[year]) / (n_prosumer)
                    else:
                        cashflow_prosumer[year] = prosumer.ec_perf_evolution['cashflow'][year] + (p2 * r3[year]) / (
                            n_prosumer)
                        cashflow[year] = cashflow[year] - (p2 * r3[year]) / (n_prosumer)

                    cashflow_cum_prosumer[year] = cashflow_cum_prosumer[year - 1] + \
                                                  cashflow_prosumer[year]
                    cashflow_cum[year] = cashflow_cum[year - 1] + cashflow[year]
                prosumer.ec_perf_evolution['cashflow_1'] = cashflow_prosumer
                prosumer.ec_perf_evolution['cashflow_cum_1'] = cashflow_cum_prosumer
                prosumer.ec_perf_evolution['NPV_1'] = npf.npv(int_rate, cashflow_prosumer)
                prosumer.ec_perf_evolution['pbp_1'] = -prosumer.ec_perf_evolution['cashflow'][0] / np.mean(
                    cashflow_prosumer[1:])

        NPV = npf.npv(int_rate, cashflow)
        pbp = outflow0 / np.mean(cashflow[1:])
        revenue_x_member = NPV / n_consumer
        revenuter_x_member_x_year = revenue_x_member / t_inv

        self.ec_perf_evolution['NPV'] = NPV
        self.ec_perf_evolution['cashflow'] = cashflow
        self.ec_perf_evolution['cashflow_cum'] = cashflow_cum
        self.ec_perf_evolution['pbp'] = pbp

        return NPV, cashflow, cashflow_cum, pbp, revenue_x_member, revenuter_x_member_x_year, r1, r2, r3, r4, c1, c2, c3, c4, c5, c6

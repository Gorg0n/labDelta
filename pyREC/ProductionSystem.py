class ProductionSystem:
    def __init__(self, id, carrier, tech,category, power, p_con, cost_kW, inc_kW, oem_cost_kW,oem_cost_kWh, inc_kWh, inc, dur_inc_kWh,
                 dur_inc, status='new', user=None, decay=0.006,
                 life_time=20, pod=None, ma=None, ts=None, gse_mode='rid',
                 aux_components=None, cogeneration=False):
        """
            to simulate distributed generation plants
            :param id: identification code
            :param carrier: main energy vector ('electricity' or 'heat')
            :param category: ('flexible' or inflex)
            :param tech: 'pv','wt','owc','biomass_system','biomass_boiler'
            :param power: nominal power [kW] #rimane da trasformare in una lista
            :param status: 'new' or 'old'
            :param user: list of ID user physically connected
            :param decay: production annual decay (for example 0.006)
            :param life_time: useful life of the plant [year]
            :param pod: point of delivery
            :param ma: electricity market area
            :param ts: transformer station
            :param gse_mode: interface mechanism with the electrical network (ssp or rid)
            :param p_con: connections point with the network
            :param cost_kW: [€/kW]
            :param oem_cost_kW: maintenance costs [€/year/kW]
            :param oem_cost_kWh: maintenance costs [€/year/kW]
            :param inc: type of incentive for example bonus50%
            :param dur_inc:duration incentive [year]
            :param aux_components: list of auxiliary components
            """
        self.id = id
        self.carrier = carrier
        self.category=category
        self.tech = tech
        self.power = power
        self.status = status
        self.user = user
        self.decay = decay
        self.life_time = life_time
        self.pod = pod
        self.ma = ma
        self.ts = ts
        self.gse_mode = gse_mode
        self.p_con = p_con

        self.cost_kW = cost_kW  # costo iniziale dell'impianto per kW include anche quello dei componenti ausiliari
        self.inc_kW = inc_kW  # incentivi che vengono calcolati sui kW installati e dati solo una volta
        self.oem_cost_kW = oem_cost_kW  # costi di manutezione per anno non include eventuale sostituzione dei componenti
        self.oem_cost_kWh=oem_cost_kWh
        self.inc_kWh = inc_kWh  # incentivi sul kWh prodotto
        self.inc = inc  # incentivi indipendenti dalla taglia e dalla produzione
        self.dur_inc = dur_inc  # durata degli incentiv indipendenti dalla taglia e dalla produzione
        self.dur_inc_kWh = dur_inc_kWh  # durata degli incentivi dipendenti dalla produzione

        self.aux_components = aux_components
        self.cogeneration=cogeneration

        self.en_perf_evolution = {}  # {'prod':production curve,'surplus': surplus production,'unmet': not satisfied demand} level 0 and level 1
        self.cost_resource = 0

        self.cost = self.power * self.cost_kW
        self.oem_cost = self.power * self.oem_cost_kW
        self.inc_kW = self.power * self.inc_kW

        if self.carrier == 'electricity':
            if self.gse_mode == 'ssp':
                if self.power <= 3:
                    gse_tax = 4 * self.p_con
                elif self.power > 3 and self.power <= 20:
                    gse_tax = 4 * self.p_con + 30
                else:
                    gse_tax = 4 * self.p_con + 30 + 1 * (self.power - 20)

            elif self.gse_mode == 'rid':
                if self.tech == 'pv':

                    if self.power <= 20:
                        gse_tax = min(0.7 * self.power, 10000)
                    elif self.power > 20 and self.power <= 200:
                        gse_tax = min(20 * 0.7 + (self.power - 20) * 0.65, 10000)
                    else:
                        gse_tax = min(
                            20 * 0.7 + (200 - 20) * 0.65 + (self.power - 200) * 0.6,
                            10000)
                elif self.tech == 'wt':

                    if self.power <= 20:
                        gse_tax = min(0.9 * self.power, 10000)
                    elif self.power > 20 and self.power <= 200:
                        gse_tax = min(20 * 0.9 + (self.power - 20) * 0.8, 10000)
                    else:
                        gse_tax = min(
                            20 * 0.9 + (200 - 20) * 0.8 + (self.power - 200) * 0.7,
                            10000)

                elif self.tech == 'owc':

                    if self.power <= 20:
                        gse_tax = min(1.2 * self.power, 10000)
                    elif self.power > 20 and self.power <= 200:
                        gse_tax = min(20 * 1.2 + (self.power - 20) * 1, 10000)
                    else:
                        gse_tax = min(
                            20 * 1.2 + (200 - 20) * 1 + (self.power - 200) * 0.9,
                            100000)
                else:
                    gse_tax = 0

            else:
                gse_tax = 0


        else:
            gse_tax = 0
            self.inc = inc

        self.gse_tax = gse_tax



import math
import numpy as np
from ProductionSystem import ProductionSystem


class WindTurbine(ProductionSystem):

    def __init__(self,id, carrier, tech='wt', power=1, status='new', user=None, decay='0.006',
                 life_time='20', pod=None, ma=None, ts=None, gse_mode='rid',
                 p_con=None, cost_kW=4000, oem_cost_kW=40, inc=None, dur_inc=None, aux_components=None, site_elv=0, site_ht=30, site_roht=46, alpha_t=0.14, rotor_ht=60, rotor_di=76,
                 sensr_ht=60, lost=0, num=1, sher_exp=0.16,
                 turb_int=0, pwr_ratd=2000, spd_ratd=15, num_pair=25, ws_data=None, p_data=None):

        super().__init__(id=id, carrier=carrier, tech=tech, power=power, status=status, user=user, decay=decay,
                 life_time=life_time, pod=pod, ma=ma, ts=ts, gse_mode=gse_mode,
                 p_con=p_con, cost_kW=cost_kW, oem_cost_kW=oem_cost_kW, inc=inc, dur_inc=dur_inc, aux_components=aux_components)
        """
        :param id: identification code
        :param carrier: energy vector 'electricity' or 'heat'
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
        :param inc: type of incentive for example bonus50%
        :param dur_inc:duration incentive [year]
        :param aux_components: list of auxiliary components
        :param site_elv: The height of the site above sea level [m].
        :param site_ht: The height above site ground level at which the wind data was collected [m].
        :param site_roht: Hub height of WECS, as installed on site [m].
        :param alpha_t: Site wind shear exponent.
        :param rotor_ht: Rotor center height [m].
        :param rotor_di: Rotor diameter [m]
        :param sensr_ht: Sensor Height for data pairs given here below (often rotor center height) [m].
        :param lost: The percentage of turbine output power that is lost due to ineffiencies and transmission.
        :param num: Number of exactly similar turbines.
        :param sher_exp: Power-law exponent for vertical wind profile.
        :param turb_int: Turbulence intensity valid for power curve.
        :param pwr_ratd: Rated power of the turbine [kW].
        :param spd_ratd: Rated wind speed [m/s].
        :param num_pair: Number of (wind speed, power) data pairs to calculate power curve.
        :param ws_data: wind speed data to calculate power curve on the manufacturer's spec sheet [m/s]
        :param p_data: measured power on the manufacturer's spec sheet [kW]
        """
        self.t_std = 288.15  # standard temperature sea-level conditions [K]
        self.b = 6.5 / 1000  # vertical thermal gradient [K/m]
        self.g = 9.8066  # gravity [m/s2]
        self.r = 287  # specific gas constant for dry air [J/(kg*K)]
        self.bp_stp = 101325  # standard pressure sea-level conditions [Pa]
        self.gamma = self.g / (self.r * self.b)
        self.air_dens = 1.225  # power curve air density [kg/m3]

        self.site_elv = site_elv
        self.site_ht = site_ht
        self.site_roht = site_roht
        self.alpha_t = alpha_t
        self.rotor_ht = rotor_ht
        self.rotor_di = rotor_di
        self.sensr_ht = sensr_ht
        self.lost = lost
        self.num = num
        self.sher_exp = sher_exp
        self.turb_int = turb_int
        self.pwr_ratd = pwr_ratd
        self.spd_ratd = spd_ratd
        self.num_pair = num_pair
        self.ws_data = ws_data
        self.p_data = p_data

        if self.alpha_t > 1:
            self.alpha_t = 1 / 7

        self.prod = None
        self.an_prod = None
        self.en_perf_evolution = {}  # {'prod':production curve,'surplus': surplus production,'no-coverage': not satisfied demand} level 0 and level 1

    def compute_power_curve(self, ws, ws_data):
        """
        :param ws:wind speed [m/s]
        :return: power [W]
        """
        j = 1
        if ws < ws_data[3]:
            p = 0
        elif ws <= ws_data[self.num_pair] and j < self.num_pair:
            while ws_data[j] <= ws and j < self.num_pair:
                j = j + 1
            p = self.p_data[j - 1] + (ws - ws_data[j - 1]) * (self.p_data[j] - self.p_data[j - 1]) / (
                    ws_data[j] - ws_data[j - 1])
        else:
            p = self.p_data[self.num_pair]

        return p

    def compute_pwecs(self, control_signal, t, ws, bp):
        """
        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power.
        :param t: ambient temperature [°C]
        :param ws: wind speed [m/s]
        :param bp:barometric pressure [Pa]
        :return: power corrected [W]
        """

        if bp < 10000:  # not use real-time data, use model
            bp = self.bp_stp * (1 - (self.b * self.site_elv / self.t_std)) ** self.gamma

        if self.sensr_ht != self.rotor_ht:
            ws_data = self.ws_data * (self.rotor_ht / self.sensr_ht) ** self.sher_exp
        else:
            ws_data = self.ws_data

        t_elv_i = t + 273.15
        rho_t = bp / (self.r * t_elv_i)
        v_hub = ws * (self.site_roht / self.site_ht) ** self.alpha_t

        p = self.compute_power_curve(ws=v_hub, ws_data=ws_data)

        pwecs = p * rho_t / self.air_dens
        newvrat = self.spd_ratd * (rho_t / self.air_dens) ** (1 / 3)

        if pwecs > self.pwr_ratd:
            pwecs = self.pwr_ratd
        elif v_hub > newvrat:
            pwecs = self.pwr_ratd

        if control_signal != 1:
            pwecs = 0


        return pwecs

    def compute_power(self, control_signal, t, ws, bp):
        """
        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power.
        :param t: ambient temperature [°C]
        :param ws: wind speed [m/s]
        :param bp: barometric pressure [Pa]
        :return: power [W]
        """
        power = np.zeros((len(t),))
        for i in range(len(power)):
            pwecs = self.compute_pwecs(control_signal=control_signal[i], t=t[i], ws=ws[i], bp=bp[i])
            if pwecs > 0:
                pwt_ideal_i = pwecs * self.num
                power[i] = pwt_ideal_i * (1 - self.lost / 100)
            else:
                power[i] = 0
        return power

    def compute_hours(self, control_signal, t, ws, bp, time):
        """
        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power.
        :param t: ambient temperature [°C]
        :param ws: wind speed [m/s]
        :param bp: barometric pressure [Pa]
        :param time: time: 1 if hourly analysis, 0.25 if quarterly analysis
        :return: hours of continuous wind turbine operation
        """
        h = np.zeros((len(t),))
        for i in range(len(h)):
            pwecs = self.compute_pwecs(control_signal=control_signal[i], t=t[i], ws=ws[i], bp=bp[i])
            h[i] = 0
            if pwecs > 0:
                h[i] = h[i] + time
            else:
                h[i] = h[i]
        return h

    def compute_cp(self, control_signal, t, ws, bp):
        """
        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power.
        :param t: ambient temperature [°C]
        :param ws: wind speed [m/s]
        :param bp: barometric pressure [Pa]
        :return cp: wind turbine’s power coefficient
        """
        cp = np.zeros((len(t),))
        for i in range(len(cp)):
            pwecs = self.compute_pwecs(control_signal=control_signal[i], t=t[i], ws=ws[i], bp=bp[i])
            t_elv_i = t[i] + 273.15
            rho_t = bp[i] / (self.r * t_elv_i)
            v_hub_i = ws[i] * (self.site_roht / self.site_ht) ** self.alpha_t
            p_den = 0.5 * rho_t * (v_hub_i ** 3)
            area = (math.pi / 4) * (self.rotor_di ** 2)
            if pwecs > 0:
                pwt_ideal = pwecs * self.num
                pwt_net = pwt_ideal * (1 - self.lost / 100)
                cp[i] = (1000 * pwt_net / self.num) / (p_den * area)
            else:
                cp[i] = 0
        return cp

    def compute_output(self, control_signal, t, ws, bp, time):
        """

        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power
        :param t:ambient temperature [°C]
        :param ws:wind speed [m/s]
        :param bp:barometric pressure [Pa]
        :param time: time: 1 if hourly analysis, 0.25 if quarterly analysis
        :return:power [W], h, cp
        """
        power = self.compute_power(control_signal=control_signal, t=t, ws=ws, bp=bp)
        h = self.compute_hours(control_signal=control_signal, t=t, ws=ws, bp=bp, time=time)
        cp = self.compute_cp(control_signal=control_signal, t=t, ws=ws, bp=bp)
        return power, h, cp

    def update_prod(self, control_signal, t, ws, bp, time):
        """
        :param control_signal: if 0 Wind Turbine is OFF, if 1 Wind Turbine is ON and providing power
        :param t:ambient temperature [°C]
        :param ws:wind speed [m/s]
        :param bp:barometric pressure [Pa]
        :param time: time: 1 if hourly analysis, 0.25 if quarterly analysis
        :return:power [W], h, cp
        """

        self.prod = self.compute_output(control_signal=control_signal, t=t, ws=ws, bp=bp,
                                        time=time)[0]
        self.an_prod = np.array(self.prod).sum()

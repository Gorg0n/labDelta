from math import log
from math import exp
from ProductionSystem import ProductionSystem
import numpy as np


class PvPanels(ProductionSystem):
    def __init__(self, id, carrier, power, p_con, cost_kW, inc_kW, oem_cost_kW,oem_cost_kWh, inc_kWh, inc, dur_inc_kWh,
                 dur_inc, tech='pv', status='new', user=None, decay=0.006,
                 life_time=20, pod=None, ma=None, ts=None, gse_mode='rid',
                 aux_components=None,category='inflex',cogeneration=False, mode_mppt=1, isc_ref=11.32, voc_ref=43.8, t_cell_ref_c=25, I_tot_ref=1000,
                 vmppt_ref=37.2, imppt_ref=10.76, mu_isc_ref=0.04,
                 mu_voc_ref=0.24, ser_cell=60, t_cell_noct_c=44, area=1.81, n_series=1, n_parallel=1,bonus50per=True):

        super().__init__(id=id, carrier=carrier, tech=tech, power=power, status=status, user=user, decay=decay,
                 life_time=life_time, pod=pod, ma=ma, ts=ts, gse_mode=gse_mode,
                 p_con=p_con, cost_kW=cost_kW, inc_kW=inc_kW,oem_cost_kW=oem_cost_kW,oem_cost_kWh=oem_cost_kWh,inc_kWh=inc_kWh, inc=inc,dur_inc_kWh=dur_inc_kWh, dur_inc=dur_inc, aux_components=aux_components,category=category,cogeneration=cogeneration)

        """
        :param id: identification code
        :param pod:point of delivery
        :param ma: electricity market area
        :param ts: transformer station
        :param tech: 'pv'
        :param user: list of user physically connected
        :param status: 'new' or 'old'
        :param power: nominal power [kW]
        :param gse_mode: interface mechanism with the electrical network (ssp or rid)
        :param decay:decay of production plants (for example 0.006)
        :param life_time:useful life of the plant [year]
        :param p_con: connections point with the network
        :param cost_kW:initial cost [€/kW]
        :param oem_cost_kW:maintenance costs [€/year/kW]
        :param inc: type of incentive for example bonus50%
        :param dur_inc: duration incentive [year]
        :param inverter_data: list of information [power_inverter kW,inverter_cost €/kW, replace_inverter(True or False),replacement_year]
        :param battery_data:list of information [capacity_battery kWh, battery_cost €/kWh,oem_battery €/kWh/year,replace_battery (True or False),replacement_year]
        :param mode_mppt: when set to 1, PV operates at its maximum power point rather than at a load voltage specified among the components inputs. When this parameter is set to 0, the PV array will operate at the load voltage input.
        :param isc_ref: [A] The module's short circuit current reported on the manufacturer's spec sheet. Reference conditions are typically 1000 W/m2 (incident solar radiation) and 25C (ambient temperature)
        :param voc_ref: [V] The module's open circuit voltage reported on the manufacturer's spec sheet.
        :param t_cell_ref_c:[°C] The module temperature at which the manufacturer reports open circuit voltage and short circuit current. According to the definition of standard test conditions, this value should be set to 25C. The definition of STC is: 1000 W/m2, module temperature 25°C, AM 1.5 after factory light soaking.
        :param I_tot_ref:[W/m^2] The solar radiation level at which the manufacturer reports open circuit voltage and short circuit current. This value is typically 1000 W/m2.
        :param vmppt_ref:[V] The module's maximum power point voltage reported on the manufacturer's spec sheet.
        :param imppt_ref:[A] The module's maximum power point current reported on the manufacturer's spec sheet.
        :param mu_isc_ref:[A/K] Temperature coefficient of Isc (ref. cond)
        :param mu_voc_ref:[V/K] Temperature coefficient of Voc (ref. cond)
        :param ser_cell: Number of cells wired in series
        :param t_cell_noct_c:[°C] Module temperature at NOCT
        :param area:[m^2] Module area
        :param n_series:Number of modules in series
        :param n_parallel:Number of modules in parallel
        :param bonus50per: tipo di incentivo
        """
        # constant
        self.qbz = 11604.45  # electron charge/Boltzmann constant [C*K/J]
        self.t_amb_noct = 293.15  # ambient temperature at NOCT [K]
        self.I_tot_noct = 800  # irradiance at NOCT [W/m^2]
        self.r_shunt = 10 ** 10  # shunt resistance
        self.ta_normal = 0.95  # transimittance-absorptance product of the cover at normal incidence
        self.eg = 1.12  # material bandgap energy (1.12 eV for Si, 1.35 eV for GaAs) [eV=J/C]

        self.mode_mppt = mode_mppt
        self.isc_ref = isc_ref
        self.voc_ref = voc_ref
        self.t_cell_ref_c = t_cell_ref_c + 273.15  # Kelvin-Celsius
        self.I_tot_ref = I_tot_ref
        self.vmppt_ref = vmppt_ref
        self.imppt_ref = imppt_ref
        self.mu_isc_ref = mu_isc_ref
        self.mu_voc_ref = mu_voc_ref
        self.ser_cell = ser_cell
        self.t_cell_noct_c = t_cell_noct_c + 273.15  # Kelvin-Celsius
        self.area = area
        self.n_series = n_series
        self.n_parallel = n_parallel
        self.bonus50per=bonus50per

        # efficienza di riferimento e coefficiente di perdite termiche dell'array
        self.eff_ref = (self.imppt_ref * self.vmppt_ref) / (self.I_tot_ref * self.area)
        self.ul = (self.I_tot_noct * self.ta_normal) / (self.t_cell_noct_c - self.t_amb_noct)  # [W/(m^2 K)]

        # compute rserie
        self.r_serie = self.compute_rserie()

        # compute gam, il_ref e io_ref from rserie
        self.gam = self.qbz * (self.vmppt_ref - self.voc_ref + self.imppt_ref * self.r_serie) / (
                self.t_cell_ref_c * log(1 - self.imppt_ref / self.isc_ref))  # [dimensionless]
        self.il_ref = self.isc_ref  # [A]
        self.io_ref = self.il_ref / exp(self.qbz * self.voc_ref / (self.gam * self.t_cell_ref_c))  # [A]

        # From single module to array
        self.il_ref = self.n_parallel * self.il_ref  # [A]
        self.io_ref = self.n_parallel * self.io_ref  # [A]
        self.gam = self.n_series * self.gam  # [adimensionale]
        self.r_serie = (self.n_series / self.n_parallel) * self.r_serie  # [ohm]
        self.a = self.gam / (self.n_series * self.ser_cell)  # [adimensionale]
        self.voc_ref = self.voc_ref * self.n_series  # [V]
        self.array_area = self.area * self.n_series * self.n_parallel  # [m^2]
        self.pmaxappx = self.vmppt_ref * self.imppt_ref * self.n_series * self.n_parallel / 1000  # [W]
        self.vmaxappx = self.vmppt_ref * self.n_series  # [V]
        self.en_perf_evolution = {}  # {'prod':production curve,'surplus': surplus production,'no-coverage': not satisfied demand} level 0 and level 1


        if self.bonus50per == True:
            self.inc += min(self.cost / 10 / 2, 96000 / 10 / 2)
        else:
            self.inc += 0

    def compute_rserie(self):
        """
        compute rserie by bisection method
        :return: rserie [ohm]
        """
        # upper limits
        rs_up = ((self.ser_cell * self.t_cell_ref_c * log(
            1 - self.imppt_ref / self.isc_ref) / self.qbz) + self.voc_ref - self.vmppt_ref) / self.imppt_ref
        gam_up = self.ser_cell
        a_up = 1
        io_up = self.isc_ref * exp(-self.qbz * self.voc_ref / (gam_up * self.t_cell_ref_c))

        # lower limits
        rs_low = 0
        gam_low = self.qbz * (self.vmppt_ref - self.voc_ref) / (
                self.t_cell_ref_c * log(1 - self.imppt_ref / self.isc_ref))
        a_low = gam_low / self.ser_cell
        io_low = self.isc_ref * exp(-self.qbz * self.voc_ref / (gam_low * self.t_cell_ref_c))

        while abs(rs_up - rs_low) > 0.0005:
            rs_new = 0.5 * (rs_up + rs_low)  # [ohm]
            gam_new = self.qbz * (self.vmppt_ref - self.voc_ref + self.imppt_ref * rs_new) / (
                    self.t_cell_ref_c * log(1 - self.imppt_ref / self.isc_ref))
            a_new = gam_new / self.ser_cell
            io_new = self.isc_ref * exp(-self.qbz * self.voc_ref / (gam_new * self.t_cell_ref_c))
            f_up = -self.mu_voc_ref + (gam_up / self.qbz) * (
                    log(1 + self.isc_ref / io_up) + (self.t_cell_ref_c / (self.isc_ref + io_up)) * (
                    self.mu_isc_ref - self.isc_ref * (
                    (self.qbz * self.eg / (a_up * self.t_cell_ref_c ** 2)) + 3 / self.t_cell_ref_c)))
            f_low = -self.mu_voc_ref + (gam_low / self.qbz) * (
                    log(1 + self.isc_ref / io_low) + (self.t_cell_ref_c / (self.isc_ref + io_low)) * (
                    self.mu_isc_ref - self.isc_ref * (
                    (self.qbz * self.eg / (a_low * self.t_cell_ref_c ** 2)) + 3 / self.t_cell_ref_c)))
            f_new = -self.mu_voc_ref + (gam_new / self.qbz) * (
                    log(1 + self.isc_ref / io_new) + (self.t_cell_ref_c / (self.isc_ref + io_new)) * (
                    self.mu_isc_ref - self.isc_ref * (
                    (self.qbz * self.eg / (a_new * self.t_cell_ref_c ** 2)) + 3 / self.t_cell_ref_c)))

            if (f_low * f_new) < 0.0:
                rs_up = rs_new

            elif (f_low * f_new) > 0.0:
                rs_low = rs_new

            gam_up = self.qbz * (self.vmppt_ref - self.voc_ref + self.imppt_ref * rs_up) / (
                    self.t_cell_ref_c * log(1 - self.imppt_ref / self.isc_ref))  # [adimensionale]
            a_up = gam_up / self.ser_cell  # [adimensionale]
            io_up = self.isc_ref * exp(-self.qbz * self.voc_ref / (gam_up * self.t_cell_ref_c))  # [A]
            gam_low = self.qbz * (self.vmppt_ref - self.voc_ref + self.imppt_ref * rs_low) / (
                    self.t_cell_ref_c * log(1 - self.imppt_ref / self.isc_ref))  # [adimensionale]
            a_low = gam_low / self.ser_cell  # [adimensionale]
            io_low = self.isc_ref * exp(-self.qbz * self.voc_ref / (gam_low * self.t_cell_ref_c))  # [A]

        r_serie = rs_new
        return r_serie

    def compute_total_radiation(self, slope, theta, I_beam, I_skydiff, I_grounddiff):
        """

        :param slope:slope of PV array [°]
        :param theta:angle of incidence of beam radiation on the array surface [°]
        :param I_beam: components of incident radiation [W/m^2]
        :param I_skydiff: components of incident radiation [W/m^2]
        :param I_grounddiff:components of incident radiation [W/m^2]
        :return:Total incident radiation [W/m^2]
        """
        if theta is not None:
            I_total=np.zeros((len(I_beam),))
            for i in range(len(I_beam)):
                I_total[i] = I_beam[i] + I_skydiff[i] + I_grounddiff[i]
                if I_total[i] > 0.1:
                    theta_diff = 59.56748 - 0.09123155 * slope - 0.00054240 * slope ** 2 + 0.00003216 * slope ** 3 - 0.00000017 * slope ** 4
                    iam_skydiffuse = 1.0 - 1.098 * 10 ** -4 * (theta_diff) + 6.26 * 10 ** -6 * (
                            theta_diff ** 2) + 6.583 * 10 ** -7 * (theta_diff ** 3) - 1.472 * 10 ** -8 * (theta_diff ** 4)
                    if iam_skydiffuse < 0.0:
                        iam_skydiffuse = 0

                    theta_gnd = 90.03182 - 0.6614549 * slope + 0.00479618 * slope ** 2 - 0.00001543 * slope ** 3 + 0.00000002 * slope ** 4
                    iam_grounddiffuse = 1.0 - 1.098 * 10 ** -4 * (theta_gnd) + 6.26 * 10 ** -6 * (
                            theta_gnd ** 2) + 6.583 * 10 ** -7 * (theta_gnd ** 3) - 1.472 * 10 ** -8 * (theta_gnd ** 4)
                    if iam_grounddiffuse < 0.0:
                        iam_grounddiffuse = 0

                    iam_beam = 1.0 - 1.098 * 10 ** -4 * (theta[i]) + 6.26 * 10 ** -6 * ((theta[i]) ** 2) + 6.583 * 10 ** -7 * (
                            (theta[i]) ** 3) - 1.472 * 10 ** -8 * ((theta[i]) ** 4)
                    if iam_beam < 0.0:
                        iam_beam = 0
                else:
                    iam_beam = 1
                    iam_grounddiffuse = 1
                    iam_skydiffuse = 1

                I_total[i] = iam_beam * I_beam[i] + iam_skydiffuse * I_skydiff[i] + iam_grounddiffuse * I_grounddiff[i]
        else:
            I_total = np.array(I_beam) + np.array(I_skydiff) + np.array(I_grounddiff)


        return (I_total)

    def compute_output_0(self, I_total, t_amb):
        """

        :param I_total:Total incident radiation [W/m^2]
        :param t_amb:ambient temperature [°C]
        :return: vmp:max. power point voltage [V]
                 imp:current at max. power [A]
                 p_max:maximum power point along IV curve [W]
                 voc:open circuit voltage [V]
                 isc:short circuit current [A]
                 t_cell:cell temperature [K]
        """
        vmp=np.zeros((len(I_total),))
        imp=np.zeros((len(I_total),))
        p_max=np.zeros((len(I_total),))
        voc=np.zeros((len(I_total),))
        isc=np.zeros((len(I_total),))
        t_cell=np.zeros((len(I_total),))
        io=np.zeros((len(I_total),))
        for i in range(len(t_amb)):
            t_amb_k_i = t_amb[i] + 273.15
            if I_total[i] < 1:
                t_cell[i] = t_amb_k_i  # Celsius-Kelvin
            else:
                t_cell[i] = t_amb_k_i + (I_total[i] * self.ta_normal - I_total[i] * self.eff_ref) / self.ul

            cellTempConv = False
            iterCount = 0
            while cellTempConv == False:
                il = (I_total[i] / self.I_tot_ref) * (
                        self.il_ref + self.mu_isc_ref * self.n_parallel * (t_cell[i] - self.t_cell_ref_c))

                if il < 0:
                    il = 0

                io[i] = self.io_ref * ((t_cell[i] / self.t_cell_ref_c) ** 3) * exp(
                    (self.qbz * self.eg / self.a) * ((1 / self.t_cell_ref_c) - (1 / t_cell[i])))

                if il > 0:
                    voc[i] = self.gam * t_cell[i] * log(il / io[i] + 1) / self.qbz
                else:
                    voc[i] = 0

                isc[i] = il

                # Newton's method for calculating the maximum power point
                if il > 0:
                    imxn = I_total[i] / self.I_tot_ref * self.n_parallel * (
                            self.imppt_ref + self.mu_isc_ref * (t_cell[i] - self.t_cell_ref_c))  # [A]
                    imxo = 0

                    while abs(imxn - imxo) > 0.001:
                        imxo = imxn
                        if imxo > il:
                            imxo = il * 0.99

                        if io[i] > 0:
                            f1 = imxo + (imxo - il - io[i]) * (log((il - imxo + io[i]) / io[i]) - imxo * self.qbz * self.r_serie / (
                                    self.gam * t_cell[i])) / (
                                         1 + (il - imxo + io[i]) * (self.qbz * self.r_serie / (self.gam * t_cell[i])))
                            f1p = 2 + (log((il - imxo + io[i]) / io[i]) - (
                                    self.qbz * self.r_serie * imxo / (self.gam * t_cell[i]))) / (
                                          (1 + (il - imxo + io[i]) * (self.qbz * self.r_serie / (self.gam * t_cell[i]))) ** 2)

                        else:
                            f1 = 0
                            f1p = 1

                        imxn = imxo - (f1 / f1p)

                    imp[i] = imxn

                    if io[i] > 0:
                        vmp[i] = log(1 + (il - imp[i]) / io[i]) * (t_cell[i] * self.gam / self.qbz) - imp[i] * self.r_serie
                    else:
                        vmp[i] = 0

                    p_max[i] = imp[i] * vmp[i]


                else:
                    imp[i] = 0
                    vmp[i] = 0
                    p_max[i] = 0

                t_cell_new = t_amb_k_i + (I_total[i] * self.ta_normal * self.array_area - p_max[i]) / (
                        self.ul * self.array_area)  # [K]

                if abs(t_cell[i] - t_cell_new) > 0.01 and iterCount < 100:
                    t_cell[i] = t_cell_new
                    iterCount = iterCount + 1
                else:
                    cellTempConv = True

        return vmp, imp, p_max, voc, isc, t_cell

    def compute_fill_factor(self, vmp, imp, voc, isc):
        """

        :param vmp:max. power point voltage [V]
        :param imp:imp:current at max. power [A]
        :param voc:open circuit voltage [V]
        :param isc:short circuit current [A]
        :return: fill factor
        """
        ff = np.zeros((len(vmp),))
        for i in range(len(vmp)):
            if voc[i] > 0 and isc[i] > 0:
                ff[i] = (vmp[i] * imp[i] / (voc[i] * isc[i]))
            else:
                ff[i] = 0
        return ff

    def compute_efficiency(self, p_max, I_total):
        """

        :param p_max:p_max:maximum power point along IV curve [W]
        :param I_total:Total incident radiation [W/m^2]
        :return: efficiency
        """
        eff = np.zeros((len(p_max),))
        for i in range(len(p_max)):
            if I_total[i] > 0:
                if self.mode_mppt > 0:
                    eff[i] = p_max[i] / (I_total[i] * self.area * self.n_series * self.n_parallel)
                else:
                    p = 0  # to implemented,missing
                    eff[i] = p / (I_total[i] * self.area * self.n_series * self.n_parallel)
            else:
                eff[i] = 0
        return eff

    def compute_output(self, slope, theta, I_beam, I_skydiff, I_grounddiff, t_amb):
        """

        :param slope: slope of PV array [°]
        :param theta:angle of incidence of beam radiation on the array surface [°]
        :param I_beam: components of incident radiation [W/m^2]
        :param I_skydiff: components of incident radiation [W/m^2]
        :param I_grounddiff:components of incident radiation [W/m^2]
        :param t_amb:ambient temperature [°C]
        :return: I_total [W/m2], vmp [V], imp [A], p_max [W], voc [V], isc [A], t_cell [°C], ff, eff
        """

        I_total = self.compute_total_radiation(slope=slope, theta=theta, I_beam=I_beam, I_skydiff=I_skydiff, I_grounddiff=I_grounddiff)
        vmp, imp, p_max, voc, isc, t_cell = self.compute_output_0(I_total=I_total, t_amb=t_amb)
        ff = self.compute_fill_factor(vmp=vmp, imp=imp, voc=voc, isc=isc)
        eff = self.compute_efficiency(p_max=p_max, I_total=I_total)


        self.en_perf_evolution['electricity'] = {}
        self.en_perf_evolution['electricity']['prod'] = p_max/1000


        return I_total, vmp, imp, p_max, voc, isc, t_cell, ff, eff



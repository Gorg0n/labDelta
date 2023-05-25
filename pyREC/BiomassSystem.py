from ProductionSystem import ProductionSystem
import numpy as np


class BiomassSystem(ProductionSystem):
    def __init__(self, id, converter, engine, cost_kW, inc_kW, oem_cost_kW, oem_cost_kWh, inc_kWh, inc, dur_inc_kWh,
                 dur_inc, category, mode, tech='BiomassSystem', status=None, decay=0.001,
                 life_time=20,
                 compressor=None, storage=None, initial_pressure_level=None,
                 gas_temperature_storage=None,
                 gas_pressure_compressor=None):

        super().__init__(id=id, carrier=engine.carrier, tech=tech, power=engine.power, status=status,
                         user=engine.user, decay=decay, category=category, cogeneration=engine.cogeneration,
                         life_time=life_time, pod=engine.pod, ma=engine.pod, ts=engine.ts, gse_mode=engine.gse_mode,
                         p_con=engine.p_con, cost_kW=cost_kW, inc_kW=inc_kW, oem_cost_kW=oem_cost_kW,
                         oem_cost_kWh=oem_cost_kWh, inc_kWh=inc_kWh,
                         inc=inc, dur_inc_kWh=dur_inc_kWh, dur_inc=dur_inc,
                         aux_components=[converter, compressor, storage, engine, converter.heat_recovery_system,
                                         engine.aux_components[0]])

        """
            :param id: identification code
            :param carrier: energy vector 'electricity' or 'heat'
            :param tech: 'pv','wt','owc','biomass_system','biomass_boiler'
            :param power: nominal power [kW]
            :param power_carrier_2: nominal power carrier 2 [kW]
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
            :param initial_pressure_level: related to storage; initial pressure/maximum allowable pressure ∈[0,1]
            :param gas_temperature_storage: temperature of gas[°C]
            :param gas_pressure_compressor: desidered output pressure [bar]
            
            """

        self.mode = mode
        self.converter = converter
        self.compressor = compressor
        self.storage = storage
        self.engine = engine
        self.mode = mode
        self.heat_recovery_system = converter.heat_recovery_system
        self.heat_recovery_system_1 = engine.aux_components[0]
        self.en_perf_evolution = {}  # {'prod':,'surplus': surplus production,'unmet': not satisfied demand}
        self.initial_pressure_level = initial_pressure_level
        self.gas_temperature_storage = gas_temperature_storage
        self.gas_pressure_compressor = gas_pressure_compressor

        aux = []
        for component in self.aux_components:
            if component != None:
                aux.append(component)
        self.aux_components = aux

    def compute_output(self, time, d=None):
        """

        :param d: electricity power demand [kW]
        :param time: 1 if hourly analysis [h], 0.25 if quarterly analysis [1/4 h]
        :return:
        """
        initial_pressure_level = self.initial_pressure_level
        gas_temperature_storage = self.gas_temperature_storage
        gas_pressure_compressor = self.gas_pressure_compressor

        if d is not None:
            len_ref = len(d)
        else:
            len_ref = len(gas_temperature_storage)

        biomass_cons = np.zeros((len_ref,))
        gas_flow_rate = np.zeros((len_ref,))
        gas_flow_rate_1 = np.zeros((len_ref,))
        power_gas = np.zeros((len_ref,))
        power_el_abs_converter = np.zeros((len_ref,))
        power_heat_developed_converter = np.zeros((len_ref,))
        power_heat_recovered_converter = np.zeros((len_ref,))
        temperature_gas_input = []
        pressure_gas_input = []
        for i in range(len_ref):
            temperature_gas_input.append(self.converter.temperature_syngas)
            pressure_gas_input.append(self.converter.pressure_syngas)

        for i in range(len_ref):
            biomass_cons[i], gas_flow_rate[i], gas_flow_rate_1[i], power_gas[i], power_el_abs_converter[i], \
                power_heat_developed_converter[i], power_heat_recovered_converter[i] = self.converter.compute_output(
                time=time)

        if self.mode == 1:
            power_el_abs_compressor, q_power_developed_compressor, temperature_gas_output = self.compressor.compute_output(
                temperature_gas_input=temperature_gas_input,
                pressure_gas_input=pressure_gas_input, gas_flow_rate_input=gas_flow_rate_1,
                pressure_gas_output=gas_pressure_compressor, time=time)
            if self.carrier=='electricity':
                n, p_set, p, surplus, unmet, power_heat_developed_engine, power_heat_recovered_engine, eta_fuel, eta_el, gas_cons, tCO2, gas_surplus, gas_deficit = self.engine.compute_output(
                    d=d + power_el_abs_compressor / 1000 + power_el_abs_converter, time=time)
            else:
                n, p_set, p, surplus, unmet, power_heat_developed_engine, power_heat_recovered_engine, eta_fuel, eta_el, gas_cons, tCO2, gas_surplus, gas_deficit = self.engine.compute_output(
                    d=d, time=time)

            pressure, pressure_level_storage, gas, gas_surplus, gas_unmet = self.storage.compute_output(
                gas_flow_rate_input=gas_flow_rate_1,
                gas_flow_rate_output=gas_cons, initial_pressure_level=initial_pressure_level,
                temperature_storage=gas_temperature_storage)
            power_el_abs = power_el_abs_converter + power_el_abs_compressor/1000

            if self.heat_recovery_system is not None and self.heat_recovery_system_1 is not None:
                power_heat_recovered = power_heat_recovered_converter + power_heat_recovered_engine
            elif self.heat_recovery_system is not None and self.heat_recovery_system_1 is None:
                power_heat_recovered = power_heat_recovered_converter
            elif self.heat_recovery_system is None and self.heat_recovery_system_1 is not None:
                power_heat_recovered = power_heat_recovered_engine
            else:
                power_heat_recovered = np.zeros((len_ref,))

            self.en_perf_evolution['electricity'] = {}
            self.en_perf_evolution['heat'] = {}
            self.en_perf_evolution['electricity']['prod'] = p
            self.en_perf_evolution['electricity']['surplus'] = surplus
            self.en_perf_evolution['electricity']['unmet'] = unmet
            self.en_perf_evolution['heat']['prod'] = power_heat_recovered

            if self.carrier == 'heat':
                self.en_perf_evolution['electricity']['prod'] = p - power_el_abs_compressor/1000 - power_el_abs_converter
                self.en_perf_evolution['heat']['surplus'] = self.engine.en_perf_evolution['heat']['surplus']
                self.en_perf_evolution['heat']['unmet'] = self.engine.en_perf_evolution['heat']['unmet']

            syngas_cons_kg = sum(gas_cons) * self.engine.fuel.rho
            biomass_cons_kg_1 = syngas_cons_kg * self.engine.fuel.lhv / (self.converter.eff * self.converter.biomass.lhv)
            biomass_cons_kg=sum(biomass_cons)
            self.cost_resource = biomass_cons_kg_1 / 1000 * self.converter.biomass.cost

            self.output = [p, surplus, unmet, power_heat_recovered, power_el_abs, power_el_abs_converter,
                           power_el_abs_compressor, n, p_set, gas_flow_rate, gas_flow_rate_1, gas_cons, gas_surplus,
                           gas_deficit, pressure, pressure_level_storage,biomass_cons,biomass_cons_kg_1]

            return p, surplus, unmet, power_heat_recovered, power_el_abs, n, p_set, gas_flow_rate, gas_flow_rate_1, gas_cons, gas_surplus, gas_deficit, pressure, pressure_level_storage

        else:
            n, p_set, p, surplus, unmet, power_heat_developed_engine, eta_fuel, eta_el, gas_cons, tCO2, gas_surplus, gas_deficit = self.engine.compute_output(
                fuel_flow_rate=gas_flow_rate, time=time)

            if self.heat_recovery_system is not None and self.heat_recovery_system_1 is not None:
                power_heat_recovered = power_heat_developed_converter + power_heat_developed_engine
            elif self.heat_recovery_system is not None and self.heat_recovery_system_1 is None:
                power_heat_recovered = power_heat_developed_converter
            elif self.heat_recovery_system is None and self.heat_recovery_system_1 is not None:
                power_heat_recovered = power_heat_developed_engine
            else:
                power_heat_recovered = np.zeros((len_ref,))

            self.en_perf_evolution['electricity'] = {}
            self.en_perf_evolution['heat'] = {}
            self.en_perf_evolution['electricity']['prod'] = p
            self.en_perf_evolution['electricity']['surplus'] = surplus
            self.en_perf_evolution['electricity']['unmet'] = unmet
            self.en_perf_evolution['heat']['prod'] = power_heat_recovered

            return p, surplus, unmet, power_heat_recovered, power_el_abs_converter, n, p_set, gas_flow_rate, gas_flow_rate_1, gas_cons, gas_surplus, gas_deficit

    def check_eff_calculate_cb(self, time, value_cb,self_cons_heat=None):

        pes, pes_car, h, cb_revenue = self.engine.check_eff_calculate_cb(time=time, value_cb=value_cb,self_cons_heat=self_cons_heat)

        return pes, pes_car, h, cb_revenue

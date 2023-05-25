# import pyREC modules
from User import User
from Prosumer import Prosumer
from Integrator import integrator
from Constants import *
from Rec import Rec
from PvPanels import PvPanels
from Fuel import Fuel
from InternalCombustionEngine import InternalCombustionEngine
from BiomassSystem import BiomassSystem
from AuxiliaryComponent import AuxiliaryComponent
from BiomassGasifier import BiomassGasifier
from HeatRecoverySystem import HeatRecoverySystem
from GasStorage import GasStorage
from GasCompressor import GasCompressor

# import other modules
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d
import numpy as np

# set simulation time step
time = 0.25  # 0.25 for quarter-hour analysis- 1 for hourly analysis

# import input data
    ## data meteo for wheather-dependent technologies
data_meteo = pd.read_csv('SdewesInput/DataMeteo/datameteo.csv', delimiter=',')
    ## characteristic curve from producer datasheet
ice_curve_1200 = pd.read_csv('SdewesInput/Ice_curve_input_1200.csv', delimiter=';')
ice_curve_600 = pd.read_csv('SdewesInput/Ice_curve_input_600.csv', delimiter=';')
ice_curve_300 = pd.read_csv('SdewesInput/Ice_curve_input_300.csv', delimiter=';')
ice_curve_200 = pd.read_csv('SdewesInput/Ice_curve_input_200.csv', delimiter=';')
ice_curve_100 = pd.read_csv('SdewesInput/Ice_curve_input_100.csv', delimiter=';')
ice_curve_50 = pd.read_csv('SdewesInput/Ice_curve_input_50.csv', delimiter=';')
ice_curve_35 = pd.read_csv('SdewesInput/Ice_curve_input_35.csv', delimiter=';')
ice_curve_19 = pd.read_csv('SdewesInput/Ice_curve_input_19.csv', delimiter=';')
    ## quarter-hour consumption curves for each vector
load_el = pd.read_csv('SdewesInput/user_qc_kW_el.csv', delimiter=';')
load_heat = pd.read_csv('SdewesInput/user_qc_kW_heat.csv', delimiter=';')
    ## physical variable
df_biomasssystem1 = pd.read_csv('TestBiomassModel/BiomassSystem_input.csv', delimiter=';')
df = pd.read_csv('TestBiomassModel/Compressor_input.csv', delimiter=';')
t_input_storage = df['t_input']


# define users
school = User(id='school', dem=[load_el['school'], load_heat['school']], pod='None', group='pubblic institue',
              plants='pv1', carrier=['electricity', 'heat'])
pmi = User(id='pmi', dem=[load_el['pmi'], load_heat['pmi']], pod='None', group='pmi', plants='pv2',
           carrier=['electricity', 'heat'])

user_potential = []
for i in range(0, 50):
    user_potential.append(
        User(id='user{0}'.format(i), dem=load_el['user{0}'.format(i)], pod='None', group='residential',
             carrier=['electricity']))


# define auxiliary components
cost_inverter_kW = 140
inverter1 = AuxiliaryComponent(id='inverter1', tech='inverter', replacement_year=10, replacement_cost=840)


# define production plants
n_series = 150
n_parallel = 1
power_pv = n_series * n_parallel * 0.4

pv1 = PvPanels(id='pv1', pod=None, ma=None, ts=None, tech='pv', user='user0', status='new', power=power_pv,
               gse_mode='rid',
               decay=0.006, life_time=20, p_con=1, cost_kW=1500, oem_cost_kW=40, inc=3750, dur_inc=10, mode_mppt=1,
               isc_ref=11.32,
               voc_ref=43.8, t_cell_ref_c=25, I_tot_ref=1000,
               vmppt_ref=37.2, imppt_ref=10.76, mu_isc_ref=0.04,
               mu_voc_ref=0.24, ser_cell=60, t_cell_noct_c=44, area=1.81, n_series=n_series, n_parallel=n_parallel,
               carrier='electricity', dur_inc_kWh=0, inc_kWh=0, inc_kW=0, aux_components=[inverter1], oem_cost_kWh=0)

power_cog = 300
syngas2 = Fuel(fuel='syngas')
biomass2 = Fuel(fuel='wood pellets')
biomassgasifier2 = BiomassGasifier(id='biomassgasifier2', biomass_flow_rate=250, en_el_abs_rate=3 / 115,
                                   heat_developed_rate=70 / 115, temperature_syngas=25, pressure_syngas=1, eff=0.93,
                                   biomass=biomass2, syngas=syngas2, replacement_year=0, replacement_cost=0)
compressor2 = GasCompressor(id='compresssor2', tech='GasCompressor', fuel=syngas2, replacement_year=0,
                            replacement_cost=0)
storage2 = GasStorage(id='gasstorage2', capacity=10, pressure_max=200, replacement_year=0, replacement_cost=0)
heat_recovery2 = HeatRecoverySystem(id='heat2', eff=1, replacement_year=0, replacement_cost=0)
ice2 = InternalCombustionEngine(id='ICE2', carrier='heat', power=power_cog, power_min=power_cog * 0.4,
                                power_max=power_cog * 1.2, n_min=1, n_max=1, fuel=syngas2,
                                datasheet_curve=ice_curve_300, category='flex', cogeneration=True, p_con=1,
                                cost_kW=0, inc_kWh=0, inc_kW=0, inc=0, dur_inc=0, dur_inc_kWh=0, oem_cost_kW=0,
                                oem_cost_kWh=0, aux_components=[heat_recovery2])
biomasssystem_mode1 = BiomassSystem(id='biomasssystem1', converter=biomassgasifier2, engine=ice2, mode=1,
                                    compressor=compressor2, storage=storage2, cost_kW=4100, inc_kW=0, oem_cost_kW=0,
                                    inc_kWh=0, inc=0, dur_inc_kWh=20, dur_inc=0, category='flex',
                                    initial_pressure_level=0.4, gas_temperature_storage=t_input_storage,
                                    gas_pressure_compressor=df_biomasssystem1['p'], oem_cost_kWh=0.025)

# calculate production plants
pv1.compute_output(slope=30, theta=None, I_beam=data_meteo['I_beam [W/m2]'], I_skydiff=data_meteo['I_skydiff [W/m2]'],
                   I_grounddiff=data_meteo['I_grounddiff [W/m2]'], t_amb=data_meteo['t_amb [C]'])


# define prosumer
prosumer1 = Prosumer(id='School complex', plant=[pv1], user=[school], list_carrier=['electricity'])

# calculate prosumer energy perfomance without REC

out_prosumer1 = prosumer1.energy_perfomance(time=0.25)
prosumer2 = Prosumer(id='pmi', plant=[biomasssystem_mode1], user=[pmi], list_carrier=['heat', 'electricity'])
out_prosumer2 = prosumer2.energy_perfomance(time=0.25)

# calculate prosumer economical perfomance without REC

pr_export_el = 232.5 * 0.9  # [€/MWh]
pr_import_el = 232.5  # [€/MWh]
pr_import_methane = 0.75  # [€/Sm3]
kwh_term_methane = 8.19  # [kWhterm] ottenuti da 1Sm3 di metano considerando un rendimento della caldaia pari a 0.91]
pr_import_heat = pr_import_methane / kwh_term_methane * conv_kWh_MWh
out_ec_prosumer2 = prosumer2.economic_perfomance(time=0.25, t_inv=20, down_payment_percentual=100, t_res=0,
                                                 int_rate=0.03,
                                                 pr_import={'electricity': pr_import_el, 'heat': pr_import_heat},
                                                 pr_export={'electricity': pr_export_el, 'heat': 0}, tax=20,
                                                 value_cb=250)

out_ec_prosumer1 = prosumer1.economic_perfomance(time=0.25, t_inv=20, down_payment_percentual=100, t_res=0,
                                                 int_rate=0.03,
                                                 pr_import={'electricity': pr_import_el},
                                               pr_export={'electricity': pr_export_el}, tax=20)

# define consumers
n_independent_var = 50
X = [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1,
     1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1]
selection = {}
for i in range(0, n_independent_var):
    selection[user_potential[i]] = X[i]
user_list = []
for user in selection:
    if selection[user] == 1:
        user_list.append(user)
for element in user_list:
    print(element.id)

# define Rec

rec = Rec(id='Rec', prosumer=[prosumer1, prosumer2], consumer=user_list, rec_plant=[], list_carrier=['electricity'])

# calculate energy perfomance REC
out_en_rec = rec.energy_perfomance(time=0.25)

#calculate economic perfomance REC
out_ec_rec = rec.economic_perfomance(time=time, t_inv=20, down_payment_percentual=100, t_res=0, int_rate=0.03,
                                     inc_shared={'electricity': 118, 'heat': 0},
                                     pr_export={'electricity': pr_export_el, 'heat': 0}, tax=20, p1=0.8, p2=0.2)


carrier_list = ['electricity']

# Calculate global energy pefomance (annual,month,day)
summary = True
if summary == True:
    annual = {}
    month = {}
    day = {}
    annual[rec.id] = {}
    month[rec.id] = {}
    for carrier in carrier_list:
        annual[rec.id][carrier] = {}
        month[rec.id][carrier] = {}
        for variable in rec.en_perf_evolution[carrier].keys():
            annual[rec.id][carrier][variable] = integrator(
                dataseries=rec.en_perf_evolution[carrier][variable], unit='power',
                period='year') / 1000
            month[rec.id][carrier][variable] = integrator(
                dataseries=rec.en_perf_evolution[carrier][variable], unit='power',
                period='month') / 1000

# Generate graphs
graph = True
if graph == True:

    # energy demand
    period_list = ['Winter', 'Spring', 'Summer', 'Autumn']
    start_list = [0, 2880, 5088, 7296]  # in hour
    resolution = 24 * 7  # in hour

    for carrier in carrier_list:
        fig3 = plt.plot()
        plt.xlabel('time [0.25 h]')
        if carrier == 'electricity':
            plt.ylabel('Power [kW]')
        else:
            plt.ylabel('Power [kWt]')
        plt.grid()
        plt.plot(rec.en_perf_evolution[carrier]['prod'], label='production', linewidth=1)
        plt.plot(rec.en_perf_evolution[carrier]['dem'], label='demand', linewidth=1)
        plt.legend()
        plt.savefig('SdewesOutput/Rec/{0}_{1}.png'.format(rec.id, carrier))
        plt.show()
        plt.clf()
        x = np.linspace(0, 8760 * int(1 / time), 8760 * int(1 / time))
        f1 = scipy.interpolate.interp1d(x, rec.en_perf_evolution[carrier]['prod'], kind='quadratic')
        f2 = scipy.interpolate.interp1d(x, rec.en_perf_evolution[carrier]['dem'], kind='quadratic')
        for start, period in enumerate(period_list):
            plt.figure()
            plt.title('{0}'.format(period))
            plt.xlabel('time [h]')
            if carrier == 'electricity':
                plt.ylabel('Power [kW]')
            else:
                plt.ylabel('Power [kWt]')
            # plt.ylim(0, 0.9 * peak * 1.2)
            plt.grid()
            produz = rec.en_perf_evolution[carrier]['prod'][
                     start_list[start] * int(1 / time):start_list[start] * int(1 / time) + resolution * int(
                         1 / time)]
            demand = rec.en_perf_evolution[carrier]['dem'][
                     start_list[start] * int(1 / time):start_list[start] * int(1 / time) + resolution * int(
                         1 / time)]
            asse = np.linspace(0, 24 * 7,
                               resolution * int(1 / time))
            plt.plot(asse, produz, label='Production',
                     color='blue', linewidth=1)
            plt.plot(asse, demand, label='Demand', color='red', linewidth=1)
            plt.fill_between(asse, produz, demand)
            if carrier == 'electricity':
                plt.fill_between(asse, produz, demand, where=(produz > demand), color='orange',
                                 interpolate=True, label='National grid export')
            else:
                plt.fill_between(asse, produz, demand, where=(produz > demand), color='orange',
                                 interpolate=True, label='Wasted')

            plt.fill_between(asse, produz, demand, where=(produz <= demand), color='lightblue',
                             interpolate=True, label='National grid import')
            plt.fill_between(asse, demand, 0, where=(produz >= demand), color='lightyellow',
                             interpolate=True)
            plt.fill_between(asse, produz, 0, where=(produz <= demand), color='lightyellow',
                             interpolate=True, label='Shared-energy')
            if period == 'Winter':
                plt.legend(loc='best')
            plt.savefig('SdewesOutput/Rec/{0}_{1}_{2}.png'.format(rec.id, carrier, period))
            plt.show()
            plt.clf()

        plt.figure(figsize=(18, 10))
        ax = plt.subplot(1, 1, 1)
        x_pos = np.arange(1, 13, 1)
        plt.title('{0}'.format(carrier), fontsize=20)
        ax.bar(x_pos - 0.24, (month[rec.id][carrier]['prod']), width=0.10, label='Production',
               color='tab:blue')
        ax.bar(x_pos - 0.12, (month[rec.id][carrier]['dem']), width=0.10, label='Demand', color='tab:orange')
        ax.bar(x_pos, (month[rec.id][carrier]['shared']), width=0.10, label='Shared energy',
               color='tab:green')
        if carrier == 'electricity':
            ax.bar(x_pos + 0.12, (month[rec.id][carrier]['surplus']), width=0.10, label='National grid export',
                   color='tab:red')
            ax.bar(x_pos + 0.24, (month[rec.id][carrier]['unmet']), width=0.10, label='National grid import',
                   color='tab:purple')
        else:
            ax.bar(x_pos + 0.12, (month[rec.id][carrier]['surplus']), width=0.10, label='Wasted',
                   color='tab:red')
            ax.bar(x_pos + 0.24, (month[rec.id][carrier]['unmet']), width=0.10, label='National grid import',
                   color='tab:purple')
        plt.legend(fontsize=20, framealpha=1, facecolor='white')
        if carrier == 'electricity':
            plt.ylabel('Energy [MWh]', fontsize=20)
        else:
            plt.ylabel('Energy [MWht]', fontsize=20)
        plt.yticks(fontsize=20)
        plt.xticks(x_pos, fontsize=20)
        plt.xlabel('Month', fontsize=20)
        plt.savefig('SdewesOutput/Rec/{0}_{1}_bar.png'.format(rec.id, carrier))
        plt.show()
        plt.clf()

out = True
# Generate output file
if out == True:
    df_total = pd.DataFrame()
    for carrier in carrier_list:
        df_carrier = pd.DataFrame(data=[annual[rec.id][carrier]], index=[carrier]).transpose()
    df_total = pd.concat([df_total, df_carrier], axis=1)

    df_ec = pd.DataFrame(data=rec.ec_perf_evolution).transpose()
    df_ec_plus = pd.DataFrame(data=[out_ec_rec[4], out_ec_rec[5]], index=['Revenuexmember', 'xmemberxyear'])
    df_ec_plus_1 = pd.DataFrame(data=[prosumer1.ec_perf_evolution['NPV_1'], prosumer1.ec_perf_evolution['pbp_1'],
                                      prosumer2.ec_perf_evolution['NPV_1'], prosumer2.ec_perf_evolution['pbp_1']],
                                index=['npv_pro1', 'pbp_por1', 'npv_pro2', 'pbp_por2'])

    df_total.to_excel('SdewesOutput/Rec/summary.xlsx')
    df_ec.to_excel('SdewesOutput/Rec/summary_ec.xlsx')
    df_ec_plus.to_excel('SdewesOutput/Rec/summary_ec_1.xlsx')
    df_ec_plus_1.to_excel('SdewesOutput/Rec/summary_ec_2.xlsx')

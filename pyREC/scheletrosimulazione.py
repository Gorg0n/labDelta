from Integrator import integrator
from User import User
from Prosumer import Prosumer
from Rec import Rec
from PvPanels import PvPanels
from AuxiliaryComponent import AuxiliaryComponent

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d

time = 1

data_meteo = pd.read_csv('Inputdelta/Datimeteo/datimeteonew.csv', delimiter=',')

load_el = pd.read_csv('Inputdelta/Copia_user_1h_kW_el.csv', delimiter=';')

school = User(id='school', dem=load_el['school'], pod='None', group='public institute',
              plants='pv1', carrier=['electricity'])

user_potential = []

for i in range(0, 5):
    user_potential.append(
        User(id='user{0}'.format(i), dem=load_el['user{0}'.format(i)], pod='None', group='residential',
             carrier=['electricity']))

inverter1 = AuxiliaryComponent(id='inverter1', tech='inverter', replacement_year=10, replacement_cost=140)

pv1 = PvPanels(id='pv1', pod=None, ma=None, ts=None, tech='pv', user='school', status='new', power=130,
               gse_mode='rid',
               decay=0.06, life_time=20, p_con=1, cost_kW=1500, oem_cost_kW=40, inc=0, dur_inc=10, mode_mppt=1,
               isc_ref=9.48,
               voc_ref=46.1, t_cell_ref_c=25, I_tot_ref=1000,
               vmppt_ref=37.0, imppt_ref=8.84, mu_isc_ref=0.04,
               mu_voc_ref=0.30, ser_cell=72, t_cell_noct_c=46, area=1.95, n_series=20, n_parallel=20,
               carrier='electricity', dur_inc_kWh=0, inc_kWh=0, inc_kW=0, aux_components=[inverter1], oem_cost_kWh=0)

pv1.compute_output(slope=15, theta=None, I_beam=data_meteo['I_beam [W/m2]'], I_skydiff=data_meteo['I_skydiff [W/m2]'],
                   I_grounddiff=data_meteo['I_grounddiff [W/m2]'], t_amb=data_meteo['t_amb [C]'])

prosumer1 = Prosumer(id='school', plant=[pv1], user=[school], list_carrier=['electricity'])

out_en_prosumer1 = prosumer1.energy_perfomance(time=time)

out_ec_prosumer1 = prosumer1.economic_perfomance(time=time, t_inv=20, down_payment_percentual=100, t_res=0,
                                                 int_rate=0.03,
                                                 pr_import={'electricity': 320},
                                                 pr_export={'electricity': 320*0.8}, tax=20)

rec = Rec(id='Rec', prosumer=[prosumer1], consumer=user_potential, rec_plant=[], list_carrier=['electricity'])

out_en_rec = rec.energy_perfomance(time=time)

out_ec_rec = rec.economic_perfomance(time=time, t_inv=20, down_payment_percentual=100, t_res=0, int_rate=0.03,
                                     inc_shared={'electricity': 118, 'heat': 0},
                                     pr_export={'electricity': 320*0.8, 'heat': 0}, tax=20, p1=0.8, p2=0.2)

list_carrier = ['electricity']
list_prosumer=[prosumer1]

summary = True
if summary == True:
    annual = {}
    month = {}
    day = {}
    annual[rec.id] = {}
    month[rec.id] = {}
    for carrier in list_carrier:
        annual[rec.id][carrier] = {}
        month[rec.id][carrier] = {}
        for variable in rec.en_perf_evolution[carrier].keys():
            annual[rec.id][carrier][variable] = integrator(
                dataseries=rec.en_perf_evolution[carrier][variable], unit='power',
                period='year') / 1000
            month[rec.id][carrier][variable] = integrator(
                dataseries=rec.en_perf_evolution[carrier][variable], unit='power',
                period='month') / 1000
    for prosumer in list_prosumer:
        annual[prosumer.id] = {}
        month[prosumer.id] = {}
        for carrier in list_carrier:
            annual[prosumer.id][carrier] = {}
            month[prosumer.id][carrier] = {}
            for variable in prosumer.en_perf_evolution[carrier].keys():
                annual[prosumer.id][carrier][variable] = integrator(
                    dataseries=prosumer.en_perf_evolution[carrier][variable], unit='power',
                    period='year') / 1000
                month[prosumer.id][carrier][variable] = integrator(
                    dataseries=prosumer.en_perf_evolution[carrier][variable], unit='power',
                    period='month') / 1000

graph = True
if graph == True:

    period_list = ['Winter', 'Spring', 'Summer', 'Autumn']
    start_list = [0, 2880, 5088, 7296]
    resolution = 24 * 7

    for carrier in list_carrier:
        fig3 = plt.plot()
        plt.title('REC')
        plt.xlabel('time [h]')
        if carrier == 'electricity':
            plt.ylabel('Power [kW]')
        else:
            plt.ylabel('Power [kWt]')
        plt.grid()
        asse = np.linspace(0, 365 * 24,
                           365*24 * int(1 / time))
        plt.plot(asse,rec.en_perf_evolution[carrier]['prod'], label='production', linewidth=1)
        plt.plot(asse,rec.en_perf_evolution[carrier]['dem'], label='demand', linewidth=1)
        plt.legend()
        plt.savefig('Tutorial1Output/Rec/{0}_{1}.png'.format(rec.id, carrier))
        plt.show()
        plt.clf()
        x = np.linspace(0, 8760 * int(1 / time), 8760 * int(1 / time))
        f1 = scipy.interpolate.interp1d(x, rec.en_perf_evolution[carrier]['prod'], kind='quadratic')
        f2 = scipy.interpolate.interp1d(x, rec.en_perf_evolution[carrier]['dem'], kind='quadratic')
        for start, period in enumerate(period_list):
            plt.figure()
            plt.title('{0}_REC'.format(period))
            plt.xlabel('time [h]')
            if carrier == 'electricity':
                plt.ylabel('Power [kW]')
            else:
                plt.ylabel('Power [kWt]')
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
            plt.legend(loc='best')
            plt.savefig('Tutorial1Output/Rec/{0}_{1}_{2}.png'.format(rec.id, carrier, period))
            plt.show()
            plt.clf()

        plt.figure(figsize=(18, 10))
        ax = plt.subplot(1, 1, 1)
        x_pos = np.arange(1, 13, 1)
        plt.title('REC', fontsize=20)
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
        plt.savefig('Tutorial1Output/Rec/{0}_{1}_bar.png'.format(rec.id, carrier))
        plt.show()
        plt.clf()



        for prosumer in list_prosumer:
            for carrier in list_carrier:
                x = np.linspace(0, 8760 * int(1 / time), 8760 * int(1 / time))
                f1 = scipy.interpolate.interp1d(x, prosumer.en_perf_evolution[carrier]['prod'], kind='quadratic')
                f2 = scipy.interpolate.interp1d(x, prosumer.en_perf_evolution[carrier]['dem'], kind='quadratic')
                for start, period in enumerate(period_list):
                    plt.figure()
                    plt.title('{0} {1}'.format(period,prosumer.id))
                    plt.xlabel('time [h]')
                    if carrier == 'electricity':
                        plt.ylabel('Power [kW]')
                    else:
                        plt.ylabel('Power [kWt]')
                    plt.grid()
                    produz = prosumer.en_perf_evolution[carrier]['prod'][
                             start_list[start] * int(1 / time):start_list[start] * int(1 / time) + resolution * int(
                                 1 / time)]
                    demand = prosumer.en_perf_evolution[carrier]['dem'][
                             start_list[start] * int(1 / time):start_list[start] * int(1 / time) + resolution * int(
                                 1 / time)]
                    asse = np.linspace(0, 24 * 7, 24 * 7 * 4)
                    plt.plot(asse, produz, label='Production',
                             color='blue', linewidth=1)
                    plt.plot(asse, demand, label='Demand', color='red', linewidth=1)
                    plt.fill_between(asse, produz, demand)
                    if carrier == 'electricity':
                        plt.fill_between(asse, produz, demand, where=(produz > demand), color='orange',
                                         interpolate=True, label='Sold to grid')
                    else:
                        plt.fill_between(asse, produz, demand, where=(produz > demand), color='orange',
                                         interpolate=True, label='Wasted')

                    plt.fill_between(asse, produz, demand, where=(produz <= demand), color='lightblue',
                                     interpolate=True, label='Bought from grid')
                    plt.fill_between(asse, demand, 0, where=(produz >= demand), color='lightyellow',
                                     interpolate=True)
                    plt.fill_between(asse, produz, 0, where=(produz <= demand), color='lightyellow',
                                     interpolate=True, label='Self-consumption')
                    if period == 'Winter':
                        plt.legend(loc='best')
                    plt.savefig('Tutorial1Output/Prosumer/{0}_{1}_{2}.png'.format(prosumer.id, carrier, period))
                    plt.show()
                    plt.clf()



                plt.figure(figsize=(18, 10))
                ax = plt.subplot(1, 1, 1)
                x_pos = np.arange(1, 13, 1)
                plt.title('{0} {1}'.format(carrier,prosumer.id), fontsize=20)
                ax.bar(x_pos - 0.24, (month[prosumer.id][carrier]['prod']), width=0.10, label='Production',
                       color='tab:blue')
                ax.bar(x_pos - 0.12, (month[prosumer.id][carrier]['dem']), width=0.10, label='Demand',
                       color='tab:orange')
                ax.bar(x_pos, (month[prosumer.id][carrier]['self_cons']), width=0.10, label='Self-consumption',
                       color='tab:green')
                if carrier == 'electricity':
                    ax.bar(x_pos + 0.12, (month[prosumer.id][carrier]['surplus']), width=0.10, label='Sold to grid',
                           color='tab:red')
                    ax.bar(x_pos + 0.24, (month[prosumer.id][carrier]['unmet']), width=0.10, label='Bought from grid',
                           color='tab:purple')
                else:
                    ax.bar(x_pos + 0.12, (month[prosumer.id][carrier]['surplus']), width=0.10, label='Wasted',
                           color='tab:red')
                    ax.bar(x_pos + 0.24, (month[prosumer.id][carrier]['unmet']), width=0.10,
                           label='National grid import',
                           color='tab:purple')
                plt.legend(fontsize=20, framealpha=1, facecolor='white')
                if carrier == 'electricity':
                    plt.ylabel('Energy [MWh]', fontsize=20)
                else:
                    plt.ylabel('Energy [MWht]', fontsize=20)
                plt.yticks(fontsize=20)
                plt.xticks(x_pos, fontsize=20)
                plt.xlabel('Month', fontsize=20)
                plt.savefig('Tutorial1Output/Prosumer/{0}_{1}_bar.png'.format(prosumer.id, carrier))
                plt.show()
                plt.clf()

                plt.figure(figsize=(18, 10))
                ax = plt.subplot(1, 1, 1)
                x_pos = np.linspace(1, 21, 21)
                x_pos1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
                plt.yticks(fontsize=20)
                plt.xticks(x_pos, x_pos1, fontsize=20)
                plt.title('Main Cost and Revenue {0}'.format(prosumer.id), fontsize=20)
                plt.ylabel(' [k€]', fontsize=20)
                ax.bar(x_pos - 0.12, out_ec_prosumer1[5] / 1000, width=0.10, label='Savings',
                       color='tab:green')
                ax.bar(x_pos, out_ec_prosumer1[4] / 1000, width=0.10, label='Revenue from sale',
                       color='tab:red')
                ax.bar(x_pos + 0.12, out_ec_prosumer1[7] / 1000, width=0.10, label='Incentives',
                       color='tab:purple')
                ax.bar(x_pos - 0.12, -out_ec_prosumer1[8] / 1000, width=0.10, label='Resource cost', color='tab:orange')
                ax.bar(x_pos, -out_ec_prosumer1[9] / 1000, width=0.10, label='O&M cost',
                       color='tab:blue')
                ax.bar(x_pos + 0.12, -out_ec_prosumer1[11] / 1000, width=0.10, label='Tax',
                       color='tab:cyan')
                # legend = plt.legend(frameon=1, fontsize=20, framealpha=0)
                plt.legend(facecolor='white', framealpha=1, fontsize=20)
                plt.xlabel('Month', fontsize=20)
                plt.savefig('Tutorial1Output/Prosumer/{0}_{1}_bar1.png'.format(prosumer.id, carrier))
                plt.show()
                plt.clf()
                plt.figure(figsize=(18, 10))
                ax = plt.subplot(1, 1, 1)
                x_pos = np.linspace(1, 21, 21)
                x_pos1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
                plt.yticks(fontsize=20)
                plt.xticks(x_pos, x_pos1, fontsize=20)
                plt.title('Cumulative Cash Flow {0}'.format(prosumer.id), fontsize=20)
                plt.ylabel(' [k€]', fontsize=20)
                ax.bar(x_pos, out_ec_prosumer1[2] / 1000, width=0.10, label='Prosumer2',
                       color='tab:blue')
                # legend = plt.legend(frameon=1, fontsize=20, framealpha=0)
                plt.legend(facecolor='white', framealpha=1, fontsize=20)
                plt.savefig('Tutorial1Output/Prosumer/{0}_{1}_bar2.png'.format(prosumer.id, carrier))
                plt.show()
                plt.clf()

out = True

if out == True:

    df_total = pd.DataFrame()
    for carrier in list_carrier:
        df_carrier = pd.DataFrame(data=[annual[rec.id][carrier]], index=[carrier]).transpose()
    df_total = pd.concat([df_total, df_carrier], axis=1)

    df_ec_rec = pd.DataFrame(data=rec.ec_perf_evolution).transpose()
    df_ec_plus = pd.DataFrame(data=[out_ec_rec[4], out_ec_rec[5]], index=['Revenuexmember', 'xmemberxyear'])
    df_total.to_excel('Tutorial1Output/Rec/summary.xlsx')
    df_ec_rec.to_excel('Tutorial1Output/Rec/summary_ec.xlsx')
    df_ec_plus.to_excel('Tutorial1Output/Rec/summary_ec_1.xlsx')

    df_total = pd.DataFrame()
    i = 0
    for prosumer in list_prosumer:
        df_prosumer = pd.DataFrame()
        df_ec_prosumer = pd.DataFrame()
        df_ec_prosumer = pd.DataFrame(data=prosumer.ec_perf_evolution).transpose()
        i += 1
        df_ec_prosumer.to_excel('Tutorial1Output/Prosumer/summary_ec_prosumer{0}.xlsx'.format(i))
        for carrier in list_carrier:
            df_carrier = pd.DataFrame(data=[annual[prosumer.id][carrier]], index=[prosumer.id]).transpose()
            df_prosumer = pd.concat([df_prosumer, df_carrier], axis=1)
        df_total = pd.concat([df_total, df_prosumer], axis=1)

    df_total.to_excel('Tutorial1Output/Prosumer/summary.xlsx')











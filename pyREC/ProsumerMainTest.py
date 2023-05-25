import pandas as pd
import pvlib
from PvPanels import PvPanels
from User import User
from Prosumer_ref import Prosumer
from Fuel import Fuel
from InternalCombustionEngine_ref import InternalCombustionEngine
#import matplotlib.pyplot as plt
from AuxiliaryComponent import AuxiliaryComponent
from BiomassBoiler import BiomassBoiler
from BiomassGasifier import BiomassGasifier
from BiomassSystem_ref import BiomassSystem
from HeatRecoverySystem import HeatRecoverySystem
from GasStorage import GasStorage
from GasCompressor import GasCompressor
from Rec import Rec

#graphical analysis as output
graph = True

#file.xlsx as output
out=True

# quarterly analysis
time = 0.25  # in hours

# import energy load in kW
data_user = pd.read_csv('RecCase1/user_qc_kW.csv', delimiter=';')

# import irradiation data from pvgis
irradiation = pvlib.iotools.pvgis.get_pvgis_hourly(41.9027835, 12.4963655, start=2019, end=2019,
                                                   raddatabase='PVGIS-SARAH2', surface_tilt=30, surface_azimuth=0,
                                                   outputformat='csv', pvcalculation=True, peakpower=0.4, loss=14,
                                                   url='https://re.jrc.ec.europa.eu/api/v5_2/')
I_beam0 = irradiation[0]['poa_direct'] #
I_skydiff0 = irradiation[0]['poa_sky_diffuse']
I_grounddiff0 = irradiation[0]['poa_ground_diffuse']
t_amb0 = irradiation[0]['temp_air']

# create the quarter-hour radiation curve
I_beam = []
I_skydiff = []
I_grounddiff = []
t_amb = []
for i in range(len(I_beam0)):
    e = I_beam0[i]
    e1 = I_skydiff0[i]
    e2 = I_grounddiff0[i]
    e3 = t_amb0[i]
    I_beam.append(e)
    I_beam.append(e)
    I_beam.append(e)
    I_beam.append(e)
    I_skydiff.append(e1)
    I_skydiff.append(e1)
    I_skydiff.append(e1)
    I_skydiff.append(e1)
    I_grounddiff.append(e2)
    I_grounddiff.append(e2)
    I_grounddiff.append(e2)
    I_grounddiff.append(e2)
    t_amb.append(e3)
    t_amb.append(e3)
    t_amb.append(e3)
    t_amb.append(e3)

# define pv plant
inverter=AuxiliaryComponent(id='inverter1', tech='inverter', replacement_year=10, replacement_cost=840)
pv1 = PvPanels(id='pv1', pod=None, ma=None, ts=None, tech='pv', user='user0', status='new', power=4, gse_mode='rid',
               decay=0.006, life_time=20, p_con=1, cost_kW=1500, oem_cost_kW=40, inc=0, dur_inc=10, mode_mppt=1, isc_ref=11.32,
               voc_ref=43.8, t_cell_ref_c=25, I_tot_ref=1000,
               vmppt_ref=37.2, imppt_ref=10.76, mu_isc_ref=0.04,
               mu_voc_ref=0.24, ser_cell=60, t_cell_noct_c=44, area=1.81, n_series=5, n_parallel=2,carrier='electricity',dur_inc_kWh=0,inc_kWh=0,inc_kW=0,aux_components=[inverter])

pv1.compute_output(slope=30, theta=None, I_beam=I_beam, I_skydiff=I_skydiff, I_grounddiff=I_grounddiff, t_amb=t_amb)

fuel1 = Fuel(fuel='methane')
df_ice=pd.read_csv('TestBiomassModel/Ice_curve_input.csv', delimiter=';')
df_ice_1=pd.read_csv('TestBiomassModel/Ice_input.csv', delimiter=';')
ice = InternalCombustionEngine(id='ICE', carrier='electricity', power=3, power_min=1, power_max=5, n_min=2,
                                n_max=5, fuel=fuel1, datasheet_curve=df_ice,category='inflex',cogeneration=True,p_con=1,cost_kW=0,inc_kW=0,inc_kWh=0,inc=0,dur_inc_kWh=0,dur_inc=0,decay=0.003,oem_cost_kW=0)

df_biomassboiler1=pd.read_csv('TestBiomassModel/BiomassBoiler_input.csv', delimiter=';')
ice.compute_output(fuel_flow_rate=df_ice_1['flow_rate'],time=0.25)

biomass=Fuel(fuel='wood pellets')
biomassboiler=BiomassBoiler(id='BiomassBoiler',carrier='heat',biomass_flow_rate_max=4.6,biomass_flow_rate_min=0.8,eff=[0.91,0.93],power_el_absorbed=[0.150,0.6],pm_emmissiosn_rate=0,co_emissions_rate=0,biomass=biomass,power=20,cost_kW=0,dur_inc=0,dur_inc_kWh=0,inc=0,inc_kWh=0,inc_kW=0,oem_cost_kW=0,category='inflex')
biomassboiler.compute_output(time=0.25,biomass_flow_rate=df_biomassboiler1['flow_rate'])

user0 = User(id='user0', dem=[data_user['user0'],data_user['user0']] , pod='IT001E61366667', group='residential', plants='pv1',carrier=['electricity','heat'])
user1 = User(id='user1', dem=[data_user['user1'],data_user['user0']] , pod='IT001E61366666', group='residential',carrier=['electricity','heat'])
user2 = User(id='user2', dem=[data_user['user2'],data_user['user0']] , pod='IT001E61379433', group='residential',carrier=['electricity','heat'])

df_biomasssystem1=pd.read_csv('TestBiomassModel/BiomassSystem_input.csv', delimiter=';')
syngas=Fuel(fuel='syngas')
biomass=Fuel(fuel='wood pellets')
biomassgasifier1=BiomassGasifier(id='biomassgasifier',biomass_flow_rate=115,en_el_abs_rate=3/115,heat_developed_rate=70/115,temperature_syngas=120,pressure_syngas=1,eff=0.93,biomass=biomass,syngas=syngas,replacement_year=0,replacement_cost=0)
compressor1=GasCompressor(id='compresssor1', tech='GasCompressor', fuel=syngas,replacement_year=0,replacement_cost=0)
storage1=GasStorage(id='GasStorage',capacity=100,pressure_max=200,replacement_year=0,replacement_cost=0)
ice1=InternalCombustionEngine(id='ICE1',carrier='electricity',power=300,power_min=120,power_max=360,n_min=0,n_max=5,fuel=syngas,datasheet_curve=df_ice,category='inflex',cogeneration=True,p_con=1,cost_kW=0,inc_kWh=0,inc_kW=0,inc=0,dur_inc=0,dur_inc_kWh=0,oem_cost_kW=0)
heat_recovery=HeatRecoverySystem(id='heat1',eff=1,replacement_year=0,replacement_cost=0)
heat_recovery1=HeatRecoverySystem(id='heat1',eff=1,replacement_year=0,replacement_cost=0)
df=pd.read_csv('TestBiomassModel/Compressor_input.csv', delimiter=';')
t_input=df['t_input']
biomasssystem_mode1=BiomassSystem(id='biomasssystem1',converter=biomassgasifier1,engine=ice1,mode=2,compressor=compressor1,storage=storage1,cost_kW=0,inc_kW=0,oem_cost_kW=0,inc_kWh=0,inc=0,dur_inc_kWh=0,dur_inc=0,category='inflex',initial_pressure_level=0.4,gas_temperature_storage=t_input,gas_pressure_compressor=df_biomasssystem1['p'])
#output=biomasssystem_mode1.compute_output(time=0.25)




# define prosumer
prosumer0 = Prosumer(id='Prosumer1', plant=[pv1], user=[user0])
out1=prosumer0.energy_perfomance(time=0.25)

out=prosumer0.economic_perfomance(time=0.25,t_inv=20,down_payment_percentual=100,t_res=0,int_rate=0.03,tax=20,pr_import={'electricity':320,'heat':150},pr_export={'electricity':288,'heat':0})
rec=Rec(id='Rec',prosumer=[prosumer0],consumer=[user0],rec_plant=[pv1])


# df1 = pd.DataFrame(data=output['electricity'])
# df2 = pd.DataFrame(data=output['heat'])
# df1.to_excel('TestProsumer/test_energy_rec' + '.xlsx')
# df2.to_excel('TestProsumer/test_energy_rec1' + '.xlsx')


output=rec.economic_perfomance(time=0.25,t_inv=20,down_payment_percentual=100,t_res=0,int_rate=0.03,tax=20,inc_shared={'electricity':320,'heat': 150},pr_export={'electricity':288,'heat': 0},p1=0.2,p2=0.3)

print(output)

#df1=pd.DataFrame(data=out1['electricity'])
#df2=pd.DataFrame(data=out1['heat_cog'])
#df3=pd.DataFrame(data=out1['heat'])
#df1.to_excel('TestProsumer/test_energy_1' + '.xlsx')
#df2.to_excel('TestProsumer/test_energy_2' + '.xlsx')
#df3.to_excel('TestProsumer/test_energy_3' + '.xlsx')

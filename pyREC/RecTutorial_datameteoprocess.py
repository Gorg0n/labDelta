import pandas as pd
import pvlib

irradiation = pvlib.iotools.pvgis.get_pvgis_hourly(41.67945, 13.57354, start=2019, end=2019,
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

df=pd.DataFrame(data=t_amb)
df.to_excel('TutorialInput/DataMeteo/t_amb'+'.xlsx')

out1=pd.DataFrame(data={'I_beam [W/m2]':I_beam,'I_skydiff [W/m2]':I_skydiff,'I_grounddiff [W/m2]':I_grounddiff,'t_amb [C]':t_amb})
out1.to_csv('TutorialInput/DataMeteo/datameteo'+'.csv')
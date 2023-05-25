import calendar
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

class LoadEstimator:

    def __init__(self, time_res, year, var_hol):
        """

        :param time_res:'hour' or 'quarterhour'
        :param year: 2022
        :param var_hol: '2022-04-18' #holidays on a variable basis
        :return:data1:{DataFrame:(8760 o 35040:5)} day hour name descr band
        :return:data2:{DataFrame:(365:5)} day hour name descr band
        """

        self.time_res = time_res
        self.year = year
        self.var_hol = var_hol

    def create_calendar(self):
        if self.time_res == 'hour':
            data = pd.date_range('{0}-01-01'.format(self.year), '{0}-01-01'.format(self.year + 1), freq='H')
            data = pd.DataFrame(data=data, columns=['data'])
            data = data.drop(data.index[8760])
            freq = 24 #number of hours per day
        elif self.time_res == 'quarterhour':
            data = pd.date_range('{0}-01-01'.format(self.year), '{0}-01-01'.format(self.year + 1), freq="0h15min")
            data = pd.DataFrame(data=data, columns=['data'])
            data = data.drop(data.index[35040])
            freq = 96 #number of quarters of hour per day
        else:
            data = 0
            freq = 0

        name = []
        for i in range(len(data)):
            my_date = data['data'][i]
            e = calendar.day_name[my_date.weekday()]
            name.append(e)

        data = pd.DataFrame(data=data, columns=['data'])
        data1 = pd.DataFrame(data=np.zeros((len(data),)), columns=['day'])
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data),)), columns=['hour'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=name, columns=['name'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data),)), columns=['descr'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data),)), columns=['band'])], axis=1)

        for i in range(len(data)):
            cella = str(data['data'][i])
            string = cella.split(" ")
            data1.iat[i, 0] = string[0]
            data1.iat[i, 1] = string[1]

        italy_holyday = ['{0}-01-01'.format(self.year), '{0}-01-06'.format(self.year), '{0}-04-25'.format(self.year),
                         '{0}-05-01'.format(self.year), '{0}-06-02'.format(self.year), '{0}-08-15'.format(self.year),
                         '{0}-11-01'.format(self.year), '{0}-12-08'.format(self.year), '{0}-12-25'.format(self.year),
                         '{0}-12-26'.format(self.year), self.var_hol,'2022-08-01','2022-08-02','2022-08-03','2022-08-04','2022-08-05','2022-08-06','2022-08-07','2022-08-08','2022-08-09','2022-08-10','2022-08-11','2022-08-12','2022-08-13','2022-08-14','2022-08-15','2022-08-16','2022-08-17','2022-08-18','2022-08-19','2022-08-20','2022-08-21','2022-08-22','2022-08-23','2022-08-24','2022-08-25','2022-08-26','2022-08-27','2022-08-28','2022-08-29','2022-08-30','2022-08-31']
        data1.loc[data1['day'].isin(italy_holyday), 'descr'] = 'H'

        for i in range(len(data1)):
            if data1['descr'][i] == 0:
                if data1['name'][i] == 'Saturday':
                    data1.iat[i, 3] = 'S'  # Saturday
                elif data1['name'][i] == 'Sunday':
                    data1.iat[i, 3] = 'H'  # Holiday
        f1 = []

        for i in range(8, 10, 1):
            if self.time_res == 'hour':
                e = '0{0}:00:00'.format(i)
                f1.append(e)
            if self.time_res == 'quarterhour':
                for m, n in enumerate(['00', '15', '30', '45']):
                    e = '0{0}:{1}:00'.format(i, n)
                    f1.append(e)

        for i in range(10, 19, 1):
            if self.time_res == 'hour':
                e = '{0}:00:00'.format(i)
                f1.append(e)
            if self.time_res == 'quarterhour':
                for m, n in enumerate(['00', '15', '30', '45']):
                    e = '{0}:{1}:00'.format(i, n)
                    f1.append(e)

        f3 = []
        for i in range(0, 7, 1):
            if self.time_res == 'hour':
                e = '0{0}:00:00'.format(i)
                f3.append(e)
            if self.time_res == 'quarterhour':
                for m, n in enumerate(['00', '15', '30', '45']):
                    e = '0{0}:{1}:00'.format(i, n)
                    f3.append(e)
        for i in range(23, 24, 1):
            if self.time_res == 'hour':
                e = '{0}:00:00'.format(i)
                f3.append(e)
            if self.time_res == 'quarterhour':
                for m, n in enumerate(['00', '15', '30', '45']):
                    e = '{0}:{1}:00'.format(i, n)
                    f3.append(e)

        data1.loc[data1['hour'].isin(f1), 'band'] = 'F1'
        data1.loc[data1['hour'].isin(f3), 'band'] = 'F3'
        aux = data1.loc[data1['hour'].isin(f1)].loc[
            data1.loc[data1['hour'].isin(f1)]['name'] == 'Saturday'].index.tolist()
        data1.iloc[aux, [4]] = 'F2'
        data1.loc[data1['descr'] == 'H', 'band'] = 'F3'
        data1.loc[data1['band'] == 0, 'band'] = 'F2'

        data2 = data1.iloc[0:len(data1):freq, :]

        return (data1, data2)

    def dailyLoad(self, x, y, peak=1):
        x = np.array(x)
        y = np.array(y)
        f = interp1d(x, peak * y,kind='cubic')
        return f

    def load_estimator_mode1(self, months_band, df_f0, df_fs, df_fh):
        """
        to estimate load curve starting from consumption f1,f2,f3 for each month of the year
        :param months_band: list of consumption f1,f2,f3 for each month of the year [kWh]
        :param df_f0: normalized usage curve (weekdays)
        :param df_fs: normalized usage curve (saturday)
        :param df_fh: normalized usage curve (sunday and holidays)
        :return: nor_load: normalize load curve
        :return: load: load curve [kWh]
        """

        data1, data2 = self.create_calendar()
        ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        if self.time_res == 'hour':
            deltat = 1  # in hours
            freq = 24
        elif self.time_res == 'quarterhour':
            deltat = 0.25
            freq = 96
        else:
            deltat = 0
            freq = 0

        endTime = int(24 / deltat)
        quartorario = np.linspace(0, 23.75, endTime)

        f0 = self.dailyLoad(df_f0['x'], df_f0['y'], df_f0['y'].max())
        fs = self.dailyLoad(df_fs['x'], df_fs['y'], df_fs['y'].max())
        fh = self.dailyLoad(df_fh['x'], df_fh['y'], df_fh['y'].max())

        normalized_load = []
        for i in range(len(data2)):
            if data2.iloc[i][3] == 'S':
                for kk in quartorario:
                    normalized_load.append(fs(kk))
            elif data2.iloc[i][3] == 'H':
                for kk in quartorario:
                    normalized_load.append(fh(kk))
            else:
                for kk in quartorario:
                    normalized_load.append(f0(kk))

        data1 = pd.concat([data1, pd.DataFrame(data=np.array(normalized_load), columns=['nor_load'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data1),)), columns=['load'])], axis=1)

        k = 0
        for m, n in enumerate(ndays):
            data_month = data1[k:n * freq + k]
            k += n * freq
            f1 = data_month.loc[data_month['band'] == 'F1']
            f2 = data_month.loc[data_month['band'] == 'F2']
            f3 = data_month.loc[data_month['band'] == 'F3']

            f1index = f1.index.tolist()
            f2index = f2.index.tolist()
            f3index = f3.index.tolist()

            x1 = months_band[m][0] / f1['nor_load'].sum()
            x2 = months_band[m][1] / f2['nor_load'].sum()
            x3 = months_band[m][2] / f3['nor_load'].sum()

            for i in f1index:
                data1.loc[i, 'load'] = data1.loc[i, 'nor_load'] * x1

            for i in f2index:
                data1.loc[i, 'load'] = data1.loc[i, 'nor_load'] * x2

            for i in f3index:
                data1.loc[i, 'load'] = data1.loc[i, 'nor_load'] * x3

        nor_load = np.array(data1['nor_load'])
        load = np.array(data1['load'])

        return (data1,nor_load, load)

    def load_estimator_mode2(self, months_consumption, df_f0, df_fs, df_fh):
        """
        to estimate load curve starting from total consumption for each month of the year
        :param months_consumption: list of total consumption for each month of the year [kWh]
        :param df_f0: normalized usage curve (weekdays)
        :param df_fs: normalized usage curve (saturday)
        :param df_fh: normalized usage curve (sunday and holidays)
        :return: nor_load: normalize load curve
        :return: load: load curve [kWh]
        """

        data1, data2 = self.create_calendar()
        ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        if self.time_res == 'hour':
            deltat = 1  # in hours
            freq = 24
        elif self.time_res == 'quarterhour':
            deltat = 0.25
            freq = 96
        else:
            deltat = 0
            freq = 0

        endTime = int(24 / deltat)
        quartorario = np.linspace(0, 23.75, endTime)

        f0 = self.dailyLoad(df_f0['x'], df_f0['y'], df_f0['y'].max())
        fs = self.dailyLoad(df_fs['x'], df_fs['y'], df_fs['y'].max())
        fh = self.dailyLoad(df_fh['x'], df_fh['y'], df_fh['y'].max())

        normalized_load = []
        for i in range(len(data2)):
            if data2.iloc[i][3] == 'S':
                for kk in quartorario:
                    normalized_load.append(fs(kk))
            elif data2.iloc[i][3] == 'H':
                for kk in quartorario:
                    normalized_load.append(fh(kk))
            else:
                for kk in quartorario:
                    normalized_load.append(f0(kk))

        data1 = pd.concat([data1, pd.DataFrame(data=np.array(normalized_load), columns=['nor_load'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data1),)), columns=['load'])], axis=1)

        k = 0
        for m, n in enumerate(ndays):
            data_month = data1[k:n * freq + k]
            k += n * freq

            monthindex = data_month.index.tolist()

            if data_month['nor_load'].sum()==0:
                x=0
            else:
                x = months_consumption[m] / data_month['nor_load'].sum()

            for i in monthindex:
                data1.loc[i, 'load'] = data1.loc[i, 'nor_load'] * x

        nor_load = np.array(data1['nor_load'])
        load = np.array(data1['load'])

        return (data1,nor_load, load)

    def load_estimator_mode3(self, annual_consumption, df_f0, df_fs, df_fh):
        """
        to estimate load curve starting from total annual consumption
        :param annual_consumption: total annual consumption [kWh]
        :param df_f0: normalized usage curve (weekdays)
        :param df_fs: normalized usage curve (saturday)
        :param df_fh: normalized usage curve (sunday and holidays)
        :return: nor_load: normalize load curve
        :return: load: load curve [kWh]
        """

        data1, data2 = self.create_calendar()
        ndays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        if self.time_res == 'hour':
            deltat = 1  # in hours
            freq = 24
        elif self.time_res == 'quarterhour':
            deltat = 0.25
            freq = 96
        else:
            deltat = 0
            freq = 0

        endTime = int(24 / deltat)
        quartorario = np.linspace(0, 23.75, endTime)

        f0 = self.dailyLoad(df_f0['x'], df_f0['y'], df_f0['y'].max())
        fs = self.dailyLoad(df_fs['x'], df_fs['y'], df_fs['y'].max())
        fh = self.dailyLoad(df_fh['x'], df_fh['y'], df_fh['y'].max())

        normalized_load = []
        for i in range(len(data2)):
            if data2.iloc[i][3] == 'S':
                for kk in quartorario:
                    normalized_load.append(fs(kk))
            elif data2.iloc[i][3] == 'H':
                for kk in quartorario:
                    normalized_load.append(fh(kk))
            else:
                for kk in quartorario:
                    normalized_load.append(f0(kk))

        data1 = pd.concat([data1, pd.DataFrame(data=np.array(normalized_load), columns=['nor_load'])], axis=1)
        data1 = pd.concat([data1, pd.DataFrame(data=np.zeros((len(data1),)), columns=['load'])], axis=1)

        x = annual_consumption / data1['nor_load'].sum()

        for i in range(len(data1)):
            data1.loc[i, 'load'] = data1.loc[i, 'nor_load'] * x

        nor_load = np.array(data1['nor_load'])
        load = np.array(data1['load'])

        return (data1,nor_load,load)

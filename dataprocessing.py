# -*- coding: utf-8 -*-
"""
Created on Wed May 20 17:27:15 2020

@author: ethan

Cleans data in the SNAPSHO/performs some preliminary analysis. 
Select rows and features to include in the output file by 
modifying "output_cols."
"""
import numpy as np
import pandas as pd
import csv
import pandasql as sql
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
import os
from sklearn.impute import SimpleImputer

PATH = 'C:/Users/ethan/Documents/MEE/Research'
os.chdir(PATH)

raw_data = pd.read_csv('snapshot_13to17.csv', encoding = "ISO-8859-1")
data = raw_data.copy(deep=True)

desired_cols = ['user_id', 'timestamp', 'label-morning Sad-Happy', 'label-morning Stressed Out-Calm Relaxed ', 'State Score',
                'duration', 'awakening_duration', 'nap_duration', 'Time_in_bed', 'study_duration', 'academic_duration', 
                'extracurricular_duration', 'call_0H-24H_mean_duration', 'call_0H-24H_mean_timestamp', 'call_0H-24H_median_duration',
                'call_0H-24H_median_timestamp', 'call_0H-24H_stdev_duration', 'call_0H-24H_stdev_timestamp', 
                'call_0H-24H_total_duration', 'call_0H-24H_total_duration_incoming', 'call_0H-24H_total_duration_outgoing',
                'screen_0H-24H_total_duration', 'screen_0H-3H_total_duration', 'screen_3H-10H_total_duration', 
                'screen_10H-17H_total_duration', 'screen_17H-24H_total_duration', 'Mobility_total_distance_a_day', 'weather_avg_cloud_cover', 
                'weather_temperature_max', 'weather_precip_probability', 'exercise_duration']

output_cols = ['user_id', 'timestamp', 'label-morning Stressed Out-Calm Relaxed ', 'duration', 'academic_duration',
                'awakening_duration', 'Time_in_bed', 'call_0H-24H_total_duration', 'screen_0H-24H_total_duration', 
                'Mobility_total_distance_a_day', 'weather_avg_cloud_cover', 'weather_temperature_max', 'weather_precip_probability',
                'exercise_duration']

label_cols = ['label-morning Sad-Happy', 'label-morning Stressed Out-Calm Relaxed ', 'State Score']

datf = data[desired_cols]

#%%
N_row, N_col = datf.shape

user_list = datf['user_id'].drop_duplicates().tolist()
user_date_list = datf[['user_id', 'timestamp']]

df_list = [datf.loc[datf['user_id'] == us_id] for us_id in user_list ]
t_arr = [len(dftt['timestamp']) for dftt in df_list ]

def get_timestamp_diffs(time_arr, date_format='%m/%d/%y'):
    date_arr = [datetime.datetime.strptime(date_str, date_format) for date_str in time_arr]
    diff_arr = [(date_arr[i+1] - date_arr[i]).days for i in range(len(date_arr) - 1)]
    return diff_arr

all_day_diffs = np.hstack([get_timestamp_diffs(e['timestamp']) for e in df_list])

print("Max Diff: ", max(all_day_diffs) )
print("Min Diff: ", min(all_day_diffs) )

def get_percent_nan(data_frame):
    nnrow, nncol = data_frame.shape
    nan_percents = [100*len(data_frame.loc[pd.isnull(data_frame[col_id]), [col_id]])/nnrow for col_id in desired_cols]
    return nan_percents

percent_nans = [100*len(datf.loc[pd.isnull(datf[col_id]), [col_id]])/N_row for col_id in desired_cols]

df_by_user_percent_nans = [get_percent_nan(df_user) for df_user in df_list ]

# Remove ID's which have insufficient amount of data 

suff_ids = user_list.copy()
# REFACTOR: If perc > 50%, replace NaN values with -1 (Replace with -50 for labels)
for n in range(len(df_by_user_percent_nans)):
    perc_arr = df_by_user_percent_nans[n]
    for perc in perc_arr:
        if(perc > 15):
            print(user_list[n])
            suff_ids.remove(user_list[n])
            break
        
        
# Uncomment below line to exclude the ID's with a large enough percent of data missing in any col
new_datf = datf.loc[datf['user_id'].isin(suff_ids)]
#new_datf = datf
#%%
def convert_timestamp_to_days(time_arr, date_format='%m/%d/%y'):
    dat0 = datetime.datetime.strptime(time_arr.to_list()[0], date_format)
    date_arr = [datetime.datetime.strptime(date_str, date_format) for date_str in time_arr]
    days_arr = [(new_date - dat0).days + 1 for new_date in date_arr]
    return days_arr

#%%
new_user_list = new_datf['user_id'].drop_duplicates().tolist()
new_df_list = [new_datf.loc[new_datf['user_id'] == us_id] for us_id in new_user_list]

new_user_list = list(range(1, len(user_list)+1))

# Simplify ID's and Timestamps
dfcnt = 0
for data_frame in new_df_list:
    nnr, nnc = data_frame.shape
    data_frame['timestamp'] = convert_timestamp_to_days(data_frame['timestamp'])
    new_id = new_user_list[dfcnt]
    data_frame['user_id'] = np.repeat(new_id, nnr)
    dfcnt += 1

# Use an imputer to fill in the missing values
imp_func = SimpleImputer(missing_values=np.nan, strategy='median')



# If perc > 50%, replace NaN values with -1 (Replace with -50 for labels)
imputed_data = []
"""
for n in range(len(df_by_user_percent_nans)):
    dfr = new_df_list[n]
    perc_arr = df_by_user_percent_nans[n]
    for p in range(len(perc_arr)):
        col_name = desired_cols[p]
        perc = perc_arr[p]
        if(perc > 25):
            if col_name in label_cols:
                dfr[col_name] = dfr[col_name].fillna(-50)
            else:
                dfr[col_name] = dfr[col_name].fillna(-1)
                
    ndfr = pd.DataFrame(imp_func.fit_transform(dfr), columns=desired_cols)
    imputed_data.append(ndfr)
"""
  
      
imputed_data = [pd.DataFrame(imp_func.fit_transform(dfr), columns=desired_cols) for dfr in new_df_list] #comment out
new_data = pd.concat(imputed_data)
#output_cols = ['user_id', 'timestamp', 'label-morning Stressed Out-Calm Relaxed ', 'duration', 'academic_duration',
#                'awakening_duration', 'Time_in_bed', 'call_0H-24H_total_duration', 'screen_0H-24H_total_duration', 
#                'Mobility_total_distance_a_day', 'weather_avg_cloud_cover', 'weather_temperature_max', 'weather_precip_probability',
#                'exercise_duration']
output_data = new_data[output_cols]

#%%
OUT_PATH = PATH + '/' + 'snap_data4.csv'
output_data.to_csv(OUT_PATH)
print("Data writtien to file: " + OUT_PATH)
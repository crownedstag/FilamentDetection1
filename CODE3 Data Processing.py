#Loading packages
import pandas as pd
import os
# csv folder path
folder_path = '/file/paper/folder'  
all_files = os.listdir(folder_path)

# Filter out .csv files
csv_files = [f for f in all_files if f.endswith('.csv')]

# Read all CSV files and merge them into a DataFrame
df_list = []
for csv_file in csv_files:
    file_path = os.path.join(folder_path, csv_file)
    df = pd.read_csv(file_path)
    df_list.append(df)

merged_df = pd.concat(df_list, axis=0, ignore_index=True)

merged_df.columns


##Load sunspot data and select the date column and sunspot number columns
sunspot = pd.read_csv('/file/paper/SN_d_tot_V2.0.csv', sep=';')
selected_columns = sunspot.iloc[:, [0, 1,2,4]]

selected_columns.columns = ['year', 'month', 'day','sunspot_number']
selected_columns['date'] = pd.to_datetime(selected_columns[['year', 'month', 'day']])


#Converting the format of a dataframe in merged_df
merged_df['date'] = pd.to_datetime(merged_df['date'])

###df_vaild, rows with index 1 and 0 
###(removing incomplete image rows with index 2)
df_vaild = merged_df[merged_df['exist'].isin([0, 1])]

###df_filaments, rows with index 1
###(removing incomplete image rows with index 2 and no filaments with index 0)
df_filaments = merged_df[merged_df['exist'].isin([1])]




#####Process df_vaild to calculate the number of filaments, their lengths,
#### and the sum of their areas for each existent filament day
grouped_exist_1 = df_vaild[df_vaild['exist'] == 1].groupby('date').agg(
    filament_count=('date', 'size'),
    total_length=('length', 'sum'),
    total_area=('area', 'sum')
).reset_index()


# Record the number of filaments, the length and the sum of the filaments for the days without filaments as 0
grouped_exist_0 = df_vaild[df_vaild['exist'] == 0].groupby('date').size().reset_index(name='filament_count')

grouped_exist_0['filament_count'] = 0
grouped_exist_0['total_length'] = 0
grouped_exist_0['total_area'] = 0

final_df = pd.concat([grouped_exist_1, grouped_exist_0], axis=0).sort_values(by='date').reset_index(drop=True)

####The processed dataset and the dataset recording the daily sunspot number were combined
merged_df_10 = pd.merge(final_df, selected_columns[['date', 'sunspot_number']], on='date', how='left')

merged_df_10.columns

merged_df_10['date'] = merged_df_10['date'].dt.date
###Finish building another dataset that will be used to explore the relationship between
#### sunspot numbers and filament characteristics
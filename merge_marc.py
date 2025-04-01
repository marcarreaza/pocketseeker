import os
import glob
import pandas as pd

def merge_csv_in_folder(folder_path):
    # Get all CSV files in the folder
    csv_files = glob.glob(os.path.join(folder_path, '*.csv'))
    # Initialize an empty list to hold DataFrames
    dataframes = []
    # Loop through each CSV file and read it into a DataFrame
    for csv_file in csv_files:
        
        df = pd.read_csv(csv_file)
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            dataframes.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
    
    if dataframes:
        # Concatenate all DataFrames into one
        merged_df = pd.concat(dataframes, axis=1)
        merged_df = merged_df.drop_duplicates(subset=['Position', 'Residue'], keep='first')
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]    

        return merged_df

def merge_csvs_from_all_folders(base_folder_path):
    # Get all subfolders in the base folder
    subfolders = [f.path for f in os.scandir(base_folder_path) if f.is_dir()]
    
    # Initialize an empty list to hold merged DataFrames
    merged_dataframes = []
    
    # Loop through each subfolder and merge CSV files
    for folder in subfolders:
        merged_df = merge_csv_in_folder(folder)
        merged_dataframes.append(merged_df)
    
    # Concatenate all merged DataFrames into one
    final_merged_df = pd.concat(merged_dataframes, ignore_index=True)
    return final_merged_df

base_folder_path = './input/'

final_df = merge_csvs_from_all_folders(base_folder_path)

# Save the final merged DataFrame to a CSV file
final_df.to_csv('merged_output.csv', index=False)

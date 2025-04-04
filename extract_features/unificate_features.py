import pandas as pd
import os

def merge_matrices(files, output_file):  
    dataframes = []

    for file in files:
        df = pd.read_csv(file)

        # Asegurar que las columnas "Res" y "File" estén presentes
        if "Res" not in df.columns or "File" not in df.columns:
            raise ValueError(f"El archivo {file} no contiene las columnas 'Res' y 'File'. Verifica los datos.")

        dataframes.append(df)

    # Realizar el merge usando "Res" y "File" como claves
    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = pd.merge(merged_df, df, on=["File", "Res"], how="outer")  # "outer" conserva todas las filas

    # Guardar el DataFrame combinado
    merged_df.to_csv(output_file, index=False)

    print(f"Las matrices se han combinado correctamente en {output_file}")

# Archivos de entrada
input_files = [
    "concativity/concavity_values.csv", 
    "core_distance/distances_to_core.csv", 
    "PSSM/pssm_and_coevolution.csv", 
    "secondary_structure/SS_features.csv"
]
output_file = "matrices_combinadas.csv"

merge_matrices(input_files, output_file)



##### Adding the Binding Site column #####
main_df = pd.read_csv(output_file)

binding_site_values = []

# Ruta donde están los archivos por proteína 
binding_data_folder = "../input"

# Para cada fila de la matriz principal
for index, row in main_df.iterrows():
    file_name = row['File']       
    residue_number = row['Res']  

    # Construir la ruta al archivo de bindings correspondiente
    path = os.path.join(binding_data_folder, file_name, "residue_matrix.csv")

    # Cargar el archivo de la proteína
    if os.path.exists(path):
        binding_df = pd.read_csv(path)

        match = binding_df[binding_df['Position'] == residue_number]

        if not match.empty:
            binding_raw = match['Binding_Site'].values[0]
        
        # Convertir 'Yes'/'No' a 1/0
        if str(binding_raw).strip().lower() == 'yes':
            binding_value = 1
        else:
            binding_value = 0 
            
    else:
        print(f"Archivo no encontrado: {path}")
        binding_value = 0  # O usar np.nan si prefieres

    binding_site_values.append(binding_value)

# Añadir la columna al dataframe principal
main_df['Binding_Site'] = binding_site_values

# Eliminar la columna "File" después de fusionar
main_df = main_df.drop(columns=["File"])

# Guardar el resultado
main_df.to_csv("matriz_con_binding.csv", index=False)
print( "Binding sites añadidos correctamente en matriz_con_binding.csv")
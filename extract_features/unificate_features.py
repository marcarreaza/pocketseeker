import pandas as pd

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
 
    # Eliminar la columna "File" después de fusionar
    merged_df = merged_df.drop(columns=["File"])

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


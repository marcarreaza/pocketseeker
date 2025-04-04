import sys
import pandas as pd
import joblib
import numpy as np

sys.path.append('../extract_features')
from concativity.concavity_feature import concavity_feature
from core_distance.core_distance_feature import core_distance
from physicochemical.physicochemical_features import physicochemical_feature
from PSSM.PSSM_feature import PSSM_feature
from SASA.sasa_features import sasa_feature
from secondary_structure.ss_feature import ss_feature
from solvent.solvent_exposure_features import solvent_feature

a = "protein.pdb"
#concavity = concavity_feature(a)
#distance_to_core = core_distance(a)
#physicochemical = physicochemical_feature (a)
#pssm = PSSM_feature (a)
#sasa = sasa_feature (a) 
#ss = ss_feature (a)
#solvent = solvent_feature (a) 

#df_final = pd.concat([
    #concavity,
    #distance_to_core["Distance_to_Core"],
    #physicochemical.drop(columns=["Position", "Residue"]),
    #pssm.drop(columns=["File", "Res"]),
    #sasa.drop(columns=["Position", "Residue"]),
    #ss.drop(columns=["Res"]),
    #solvent.drop(columns=["Position", "Residue"]),
#], axis=1)

#df_final = df_final.drop(columns=["File"])
#df_final.to_csv("output.csv", index=False)
df_final = pd.read_csv('output.csv') 

## Aplicación del modelo
model = joblib.load('../random_forest_binding_site_model.joblib')
X = df_final.drop(columns=['Res'])

X.replace('-', np.nan, inplace=True)
X = pd.get_dummies(X, columns=['SS'])

for col in model.feature_names_in_:
    if col not in X.columns:
        X[col] = 0

# Reordenar columnas para que coincidan exactamente
X = X[model.feature_names_in_]

y_pred = model.predict(X)

df_predicted = {
    'Res' : df_final["Res"],
    'Binding_sites' : y_pred
}

df_predicted.to_csv("binding_sites_predictions.csv", index=False)
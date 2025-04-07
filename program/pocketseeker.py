import sys
import pandas as pd
import joblib
import numpy as np

sys.path.append('../extract_features')
from model_features import extract_features

file = sys.argv[1]
df_final = extract_features(file)

## Aplicación del modelo
model = joblib.load('../model/random_forest_binding_site_model.joblib')
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

df_predicted = pd.DataFrame(df_predicted)
df_predicted.to_csv("binding_sites_predictions.csv", index=False)
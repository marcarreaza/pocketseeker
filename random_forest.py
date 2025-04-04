import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from joblib import dump


# 1. Preparación de los datos
df = pd.read_csv("extract_features/matriz_con_binding.csv")

# Separar features (X) y etiquetas (y)
X = df.drop(columns=['Res', 'Binding_Site']) 
y = df['Binding_Site']

# Reemplazar '-' por NaN
X.replace('-', np.nan, inplace=True)

# Convierte la columna 'SS' a numérica
X = pd.get_dummies(X, columns=['SS'])


# 2. División en train y test
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# 3. Entrenamiento del modelo
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)


# 4. Predicciones
y_pred = clf.predict(X_test)

# 5. Métricas
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy del modelo: {accuracy:.4f}")

print("Matriz de confusión:")
print(confusion_matrix(y_test, y_pred))

print("Reporte de clasificación:")
print(classification_report(y_test, y_pred))

# Guardar el modelo entrenado
dump(clf, 'random_forest_binding_site_model.joblib')

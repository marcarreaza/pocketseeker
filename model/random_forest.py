import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from joblib import dump
import os


def random_forest(df, output_folder):
    df = df.dropna(subset=['Binding_Site'])

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
    clf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
    clf.fit(X_train, y_train)


    # 4. Predicciones
    #y_pred = clf.predict(X_test)
    y_probs = clf.predict_proba(X_test)[:, 1]  # probabilidades para la clase 1
    y_pred_adjusted = (y_probs > 0.3).astype(int)


    # Metrics
    accuracy = accuracy_score(y_test, y_pred_adjusted)
    conf_matrix = confusion_matrix(y_test, y_pred_adjusted)
    class_report = classification_report(y_test, y_pred_adjusted)

    # File path for storing the metrics
    metrics_path = os.path.join(output_folder, "metrics.txt")

    # Save metrics to the file
    with open(metrics_path, 'w') as f:
        f.write(f"Accuracy of the model: {accuracy:.4f}\n\n")
        f.write("Confusion Matrix:\n")
        f.write(f"{conf_matrix}\n\n")
        f.write("Classification Report:\n")
        f.write(f"{class_report}\n")

    print(f"Metrics saved to {metrics_path}")

    # Save the trained model
    model_path = os.path.join(output_folder, "random_forest_binding_site_model.joblib")
    dump(clf, model_path)
    print(f"Model saved to {model_path}")
    
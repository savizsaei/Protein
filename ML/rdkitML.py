from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.optimizers import Adam

# Step 1: Prepare molecular data (example SMILES + labels for binary classification)
smiles = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",    # aspirin - label 1 (active)
    "Cn1cnc2n(C)c(=O)n(C)c(=O)c12", # caffeine - label 0 (inactive)
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" # ibuprofen - label 1 (active)
]
labels = [1, 0, 1]

# Step 2: Generate molecular fingerprints as feature vectors
def mol_to_fp(mol):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    arr = np.zeros((1,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

mols = [Chem.MolFromSmiles(smi) for smi in smiles]
features = np.array([mol_to_fp(m) for m in mols])

# Step 3: Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.33, random_state=42)

# --- Traditional ML: Random Forest Classifier ---
rf_clf = RandomForestClassifier(n_estimators=100, random_state=42)
rf_clf.fit(X_train, y_train)
y_pred_rf = rf_clf.predict(X_test)
print("--- Random Forest Classification Report ---")
print(classification_report(y_test, y_pred_rf))


# --- Deep Learning: Feedforward Neural Network with TensorFlow/Keras ---
model = Sequential([
    Dense(512, input_shape=(2048,), activation="relu"),
    Dropout(0.3),
    Dense(256, activation="relu"),
    Dropout(0.3),
    Dense(1, activation="sigmoid")  # binary classification
])

model.compile(
    optimizer=Adam(learning_rate=0.001),
    loss="binary_crossentropy",
    metrics=["accuracy"]
)

model.fit(X_train, np.array(y_train), epochs=30, batch_size=2, verbose=0)

# Evaluate
loss, accuracy = model.evaluate(X_test, np.array(y_test), verbose=0)
print(f"--- Neural Network Test Accuracy: {accuracy:.2f} ---")

# Predict classes
y_pred_dl = (model.predict(X_test) > 0.5).astype(int)
print("--- Neural Network Classification Report ---")
print(classification_report(y_test, y_pred_dl))


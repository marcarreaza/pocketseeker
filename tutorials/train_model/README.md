## Train Your Own Model

In this tutorial, you’ll train PocketSeeker using your own protein structures and annotated binding sites.

### 🛠️ Train the model

From within the main folder, to train the model use the following command:

```bash
python3 pocketseeker.py -train_model -i tutorial/train_model/input -i tutorial/train_model/results
```

### 📁 Contents of the input data

Each subfolder represents one protein and should include:

- `protein.pdb` or `protein.mol2`: Full 3D structure of the protein.

And **one** of the following options for binding site annotation:

- `site.pdb` or `site.mol2`: A file containing only the binding site residues.
  
  **OR**

- `residues.txt`: A text file with binding residues in the following format:

  ```
  1	A
  4	V
  6	L
  ```

> 📌 Only **PDB** and **MOL2** formats are supported.

> You can mix `site.pdb` and `residues.txt` formats across different proteins in your training dataset.



### 📂 Output

By default, a results/ folder will be created in your working directory. This folder will contain:
- **`random_forest_binding_site_model.joblib`**  
  The trained Random Forest model, ready to make predictions.

- **`metrics.txt`**  
  A summary of the model’s performance, including accuracy and other evaluation metrics.

- **`features/`** *(subfolder)*  
  Contains:
  - **`total_features.csv`**  
    A combined table of features extracted from all processed proteins.
  - **One CSV file per protein**  
    Each file contains the detailed features extracted from that specific protein.

---

Once trained, you can use the model with **PocketSeeker** to predict binding sites in new protein structures.

---

Enjoy! 🧬✨


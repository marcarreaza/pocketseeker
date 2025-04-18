## Train Your Own Model

In this tutorial, you’ll train PocketSeeker using your own protein structures and annotated binding sites.

### 📁 Contents of `train_your_own/`

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

### 🛠️ Train the model

From within the `train_your_own/` folder, use the following command:

```bash
python pocketseeker.py -train_model /input
```

### 📂 Output

- A trained model will be created and saved in the working directory.
- You may see logs or summary files depending on your implementation (e.g., accuracy reports, model files).

Once trained, the model can be used with PocketSeeker on new protein inputs just like in the prediction tutorial.

---

Enjoy exploring protein binding sites with **PocketSeeker**! 🧬✨

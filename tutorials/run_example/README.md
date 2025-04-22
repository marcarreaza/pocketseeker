## Run Example Protein Prediction

In this tutorial, you'll predict the ligand-binding sites of a sample protein structure.

**Contents of the folder run_example/**
- `1a2n_1.pdb`: 3D structure of the 1a2n_1 protein used for testing.

Use the following command to run PocketSeeker on the example protein:

```bash
python3 pocketseeker.py tutorial/run_example/1a2n_1.pdb -o tutorial/run_example/results
```

### 📂 Output

By default, a results/ folder will be created in your working directory. This folder will contain:
- binding_sites_predictions.csv: The residue binding site predicted.
- modified_1a2n_1.pdb: A pdb file que se ha cambiado el B-faxctor para que sea compatible a que se pueda visualizar en el chimera.


Example output (`results/binding_sites_predictions.csv`):
```
Res	Binding_sites	Score
1	No	            0.03
2	No	            0.04
3	No	            0.07
4	No	            0.05
```

Each number corresponds to the position of a predicted ligand-binding residue.

## Run Example Protein Prediction

In this tutorial, you'll predict the ligand-binding sites of a sample protein structure.

**Contents of the folder run_example/**
- `a7k_4.pdb`: 3D structure of the a7k_4 protein used for testing.

Use the following command to run PocketSeeker on the example protein:

```bash
python pocketseeker.py /tutorial/run_example/protein.pdb
```

- A folder named `results/` will be created automatically in the current working directory.
- Inside the `results/` folder, a file named `a7k_4.txt` will be generated.
- This file contains the list of predicted binding site residues for the input protein, one per line.

Example output (`results/a7k_4.txt`):

```
45
78
102
150
```

Each number corresponds to the position of a predicted ligand-binding residue.

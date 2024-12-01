import os
import numpy as np
from rdkit.Chem import AllChem as Chem
import argparse


def types_from_file(protein_list, output_file, ligand_dir, gninatypes_dir, barycenter_dir):
    """
    Generates .types file for input to deep learning models.

    Args:
        protein_list (list): List of protein identifiers from train/test files.
        output_file (file): Output .types file to write data.
        ligand_dir (str): Directory containing ligand .sdf files.
        gninatypes_dir (str): Directory containing .gninatypes files.
        barycenter_dir (str): Directory containing barycenter .txt files.
    """
    distance = 4  # Distance threshold for labeling
    for line in protein_list:
        prot = line.strip()  # Get protein name without newline characters

        # Paths to individual files
        ligand_file = os.path.join(ligand_dir, f"{prot}_ligand.sdf")
        gninatypes_file = os.path.join(gninatypes_dir, f"{prot}_protein_nowat.gninatypes")
        barycenter_file = os.path.join(barycenter_dir, f"{prot}_bary_centers.txt")

        # Check if files exist
        if not os.path.exists(ligand_file):
            print(f"Ligand file not found: {ligand_file}")
            continue
        if not os.path.exists(gninatypes_file):
            print(f"Gninatypes file not found: {gninatypes_file}")
            continue
        if not os.path.exists(barycenter_file):
            print(f"Barycenter file not found: {barycenter_file}")
            continue

        # Load ligand
        mol = Chem.MolFromMolFile(ligand_file, sanitize=False)
        if mol is None:
            print(f"Could not load ligand: {ligand_file}")
            continue
        c = mol.GetConformer()
        atom_positions = c.GetPositions()

        # Load barycenters
        centers = np.loadtxt(barycenter_file, comments="#")
        if centers.ndim == 1:
            centers = np.expand_dims(centers, axis=0)
        if centers.shape[0] == 0:
            print(f"No barycenters found for: {barycenter_file}")
            continue

        # Process each barycenter
        for center in centers:
            label = 0
            distances = np.linalg.norm(atom_positions - center, axis=1)
            if np.any(distances <= distance):
                label = 1

            # Write to .types file
            output_file.write(
                f"{label} {center[0]} {center[1]} {center[2]} {gninatypes_file}\n"
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate .types files for deep learning models.")
    parser.add_argument('--train_list', required=True, help="Path to train.txt file.")
    parser.add_argument('--test_list', required=True, help="Path to test.txt file.")
    parser.add_argument('--ligand_dir', required=True, help="Directory containing ligand .sdf files.")
    parser.add_argument('--gninatypes_dir', required=True, help="Directory containing .gninatypes files.")
    parser.add_argument('--barycenter_dir', required=True, help="Directory containing barycenter .txt files.")
    parser.add_argument('--train_types', required=True, help="Output path for train.types file.")
    parser.add_argument('--test_types', required=True, help="Output path for test.types file.")
    
    args = parser.parse_args()

    # Read train and test protein names
    with open(args.train_list, 'r') as train_file:
        train_prots = train_file.readlines()
    with open(args.test_list, 'r') as test_file:
        test_prots = test_file.readlines()

    # Generate train.types
    with open(args.train_types, 'w') as holo_train:
        types_from_file(train_prots, holo_train, args.ligand_dir, args.gninatypes_dir, args.barycenter_dir)

    # Generate test.types
    with open(args.test_types, 'w') as holo_test:
        types_from_file(test_prots, holo_test, args.ligand_dir, args.gninatypes_dir, args.barycenter_dir)

    print("Train and test .types files generated successfully.")

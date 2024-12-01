import sys
import os
import struct
import numpy as np
sys.path.append('/usr/local/miniconda/lib/python3.12/site-packages')
import molgrid


def gninatype(file):
    """
    Creates a GNINATYPES file for model input.
    Args:
        file (str): Path to the protein PDB file.
    Returns:
        str: Path to the created GNINATYPES file.
    """
    # Create TYPES file for model input
    f = open(file.replace('.pdb', '.types'), 'w')
    f.write(file)
    f.close()

    # Replace __file__ with the current working directory
    atom_map = molgrid.FileMappedGninaTyper('/content/drive/MyDrive/MP/gninamap.txt')
    dataloader = molgrid.ExampleProvider(atom_map, shuffle=False, default_batch_size=1)
    train_types = file.replace('.pdb', '.types')
    dataloader.populate(train_types)

    example = dataloader.next()
    coords = np.array(example.coord_sets[0].coords)
    types = np.array(example.coord_sets[0].type_index)
    types = np.int_(types)

    # Write GNINATYPES file
    gninatypes_file = file.replace('.pdb', '.gninatypes')
    with open(gninatypes_file, 'wb') as fout:
        for i in range(coords.shape[0]):
            fout.write(struct.pack('fffi', coords[i][0], coords[i][1], coords[i][2], types[i]))

    os.remove(train_types)
    return gninatypes_file


def create_types(barycenter_file, protein_file, output_dir):
    """
    Creates a TYPES file for model predictions based on the barycenter and protein file.
    Args:
        barycenter_file (str): Path to the barycenter file.
        protein_file (str): Path to the protein GNINATYPES file.
        output_dir (str): Directory to save the TYPES file.
    Returns:
        str: Path to the created TYPES file.
    """
    # Check if barycenter_file exists
    if not os.path.isfile(barycenter_file):
        print(f"Error: {barycenter_file} is not a valid file.")
        return None

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create TYPES file for model predictions
    types_file = os.path.join(output_dir, os.path.basename(barycenter_file).replace('.txt', '.types'))
    with open(barycenter_file, 'r') as fin, open(types_file, 'w') as fout:
        for line in fin:
            fout.write(' '.join(line.split()) + ' ' + protein_file + '\n')

    return types_file


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <protein_folder> <barycenters_folder> <types_folder>")
        sys.exit(1)

    # Input paths
    protein_folder = sys.argv[1]  # Folder containing cleaned PDB files
    barycenters_folder = sys.argv[2]  # Folder containing barycenter files
    types_folder = sys.argv[3]  # Folder to save TYPES files

    # Ensure the barycenters directory exists
    if not os.path.isdir(barycenters_folder):
        print(f"Error: {barycenters_folder} is not a valid directory.")
        sys.exit(1)

    # Ensure the types folder exists
    if not os.path.exists(types_folder):
        os.makedirs(types_folder)

    # Loop through each cleaned PDB file in the protein folder
    for protein_file in os.listdir(protein_folder):
        if protein_file.endswith('.pdb'):
            protein_path = os.path.join(protein_folder, protein_file)

            # Generate GNINATYPES file
            gninatypes_file = gninatype(protein_path)
            print(f"GNINATYPES file created: {gninatypes_file}")

            # Process each barycenter file in the barycenters folder
            for barycenter_file in os.listdir(barycenters_folder):
                if barycenter_file.endswith('_bary_centers.txt'):
                    barycenter_path = os.path.join(barycenters_folder, barycenter_file)

                    # Generate TYPES file and save in the specified folder
                    types_file = create_types(barycenter_path, gninatypes_file, types_folder)
                    if types_file:
                        print(f"TYPES file created: {types_file}")
                    else:
                        print(f"Failed to create TYPES file for {barycenter_path}")

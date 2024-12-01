from Bio.PDB import PDBParser, PDBIO, Select
import Bio
import os
import sys

class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if Bio.PDB.Polypeptide.is_aa(residue, standard=True) else 0

def clean_pdb(input_file, output_file):
    pdb = PDBParser().get_structure("protein", input_file)
    io = PDBIO()
    io.set_structure(pdb)
    io.save(output_file, NonHetSelect())

def clean_pdb_in_folder(input_folder, output_folder):
    # Check if output folder exists, if not create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Loop through all files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".pdb"):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, f"cleaned_{filename}")
            clean_pdb(input_file, output_file)
            print(f"Processed: {filename}")

if __name__ == '__main__':
    input_folder = sys.argv[1]  # First argument is the input folder
    output_folder = sys.argv[2]  # Second argument is the output folder
    clean_pdb_in_folder(input_folder, output_folder)

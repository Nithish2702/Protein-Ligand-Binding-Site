import os
import numpy as np
import sys
import re

def get_centers(input_dir, output_dir):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bary = open(os.path.join(output_dir, 'bary_centers.txt'), 'w')
    
    # Iterate over files in the input directory
    for d in os.listdir(input_dir):
        centers = []
        masses = []
        
        if d.endswith('vert.pqr'):
            num = int(re.search(r'\d+', d).group())
            with open(os.path.join(input_dir, d)) as f:
                for line in f:
                    if line.startswith('ATOM'):
                        # Extract coordinates and mass
                        center = list(map(float, re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", ' '.join(line.split()[5:]))))[:3]
                        mass = float(line.split()[-1])
                        centers.append(center)
                        masses.append(mass)
            
            # Convert lists to numpy arrays for processing
            centers = np.asarray(centers)
            masses = np.asarray(masses)

            # Calculate the center of mass (barycenter)
            xyzm = (centers.T * masses).T
            xyzm_sum = xyzm.sum(axis=0)
            cg = xyzm_sum / masses.sum()

            # Write the results to the output file
            bary.write(f"{num}\t{cg[0]}\t{cg[1]}\t{cg[2]}\n")
    
    bary.close()

if __name__ == '__main__':
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    get_centers(input_folder, output_folder)

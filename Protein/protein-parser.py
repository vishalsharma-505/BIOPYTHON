import os
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.SeqUtils import seq1
import warnings
from Bio import BiopythonWarning

# Ignore Biopython structural warnings for cleaner terminal output
warnings.simplefilter('ignore', BiopythonWarning)

def analyze_protein(pdb_id, download_dir="PDB_Files"):
    """
    Downloads an mmCIF file, parses its SMCRA architecture, extracts the sequence,
    calculates the Center of Mass, and analyzes structural flexibility.
    """
    # 1. Directory Setup & Download
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
        
    print(f"\n[*] Fetching PDB ID: {pdb_id.upper()}...")
    pdb_list = PDBList()
    cif_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format="mmCif")
    
    # 2. Parse Structure
    print(f"[*] Parsing structure from {cif_file}...")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id.upper(), cif_file)
    
    # Variables for Analytics
    total_atoms, total_residues = 0, 0
    chain_info, b_factors, coordinates, masses = [], [], [], []
    
    print("\n" + "="*30)
    print("      STRUCTURE SUMMARY      ")
    print("="*30)
    
    # Assuming first model (standard for most X-ray/Cryo-EM structures)
    model = structure[0] 
    
    # 3. Traverse SMCRA Architecture
    for chain in model:
        chain_id = chain.id
        residues = list(chain.get_residues())
        total_residues += len(residues)
        
        # Sequence Extraction
        seq = ""
        for res in residues:
            # Filter out water molecules and ligands (heteroatoms start with 'W' or 'H')
            if res.id[0] == ' ':
                res_name = res.get_resname()
                seq += seq1(res_name) # Converts 3-letter (ALA) to 1-letter (A)
        
        chain_info.append((chain_id, len(residues), seq))
        
        # Atom level data collection
        for residue in residues:
            for atom in residue:
                total_atoms += 1
                b_factors.append((atom.get_bfactor(), atom.get_full_id(), atom.get_name()))
                coordinates.append(atom.get_coord())
                masses.append(atom.mass)

    # Print Summary
    print(f"Total Chains   : {len(model)}")
    print(f"Total Residues : {total_residues}")
    print(f"Total Atoms    : {total_atoms}")
    for cid, clen, _ in chain_info:
        print(f" -> Chain {cid}: {clen} residues")
        
    print("\n" + "="*30)
    print("       ADVANCED METRICS      ")
    print("="*30)
    
    # 4. Calculate Center of Mass (Mass-weighted average of coordinates)
    coord_matrix = np.array(coordinates)
    mass_matrix = np.array(masses)
    center_of_mass = np.average(coord_matrix, axis=0, weights=mass_matrix)
    print(f"Center of Mass (X, Y, Z): {np.round(center_of_mass, 2)} Å")
    
    # 5. Analyze B-Factors (Higher B-factor = More flexible/disordered region)
    b_factors.sort(reverse=True, key=lambda x: x[0])
    most_flexible = b_factors[0]
    # format of full_id: (id, model, chain, residue, atom)
    chain_name = most_flexible[1][2]
    res_num = most_flexible[1][3][1]
    print(f"Highest Flexibility     : Atom {most_flexible[2]} in Chain {chain_name}, Residue {res_num}")
    print(f"Max B-Factor Value      : {np.round(most_flexible[0], 2)}")
    
    # 6. Save Extracted Sequence to FASTA
    fasta_filename = f"{pdb_id.upper()}_extracted.fasta"
    with open(fasta_filename, "w") as f:
        for cid, _, cseq in chain_info:
            if cseq: # Ignore empty sequences (like purely ligand chains)
                f.write(f">{pdb_id.upper()}:Chain_{cid}\n{cseq}\n")
    print(f"\n[*] Chain sequences saved to: {fasta_filename}")

if __name__ == "__main__":
    # You can change this to any PDB ID (e.g., "1CRN", "4HHB")
    target_pdb = input("Enter a 4-letter PDB ID (e.g., 1AI0): ").strip()
    analyze_protein(target_pdb)
import os
import warnings
from Bio import BiopythonWarning
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.SeqUtils import seq1
from Bio import Align

# Ignore structural warnings for clean terminal output
warnings.simplefilter('ignore', BiopythonWarning)

def extract_sequence_from_mmcif(pdb_id, download_dir="PDB_Files"):

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
        
    print(f"[*] Fetching structure for {pdb_id.upper()}...")
    pdb_list = PDBList()
    cif_file = pdb_list.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format="mmCif")
    
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id.upper(), cif_file)
    
    # We take the first model and its first chain
    model = structure[0]
    first_chain = list(model.get_chains())[0]
    
    seq = ""
    for res in first_chain:
        # Ignore heteroatoms (water, ligands)
        if res.id[0] == ' ':
            res_name = res.get_resname()
            seq += seq1(res_name) # Convert to 1-letter amino acid code
            
    print(f"[*] Extracted sequence for {pdb_id.upper()} (Chain {first_chain.id}): Length = {len(seq)}")
    return seq

def align_sequences(seq1, seq2):
    """
    Performs pairwise alignment using custom match/mismatch/gap scores.
    """
    print("\n[*] Initializing PairwiseAligner (Smith-Waterman style)...")
    aligner = Align.PairwiseAligner()
    
    aligner.mode = 'local'
    
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    
    print("[*] Performing sequence alignment...")
    alignments = aligner.align(seq1, seq2)
    
    # Get the best alignment (the first one)
    best_alignment = alignments[0]
    
    print("\n=== Alignment Result ===")
    print(f"Alignment Score: {best_alignment.score}\n")
    print(best_alignment)

if __name__ == "__main__":
    print("--- 3D Structure to 1D Sequence Alignment Tool ---")
    
    pdb_a = input("Enter first PDB ID (e.g., 4HHB): ").strip()
    pdb_b = input("Enter second PDB ID (e.g., 1IGY): ").strip()
    
    try:
        # Step 1: Extract sequences
        sequence_a = extract_sequence_from_mmcif(pdb_a)
        sequence_b = extract_sequence_from_mmcif(pdb_b)
        
        # Step 2: Align sequences
        if sequence_a and sequence_b:
            align_sequences(sequence_a, sequence_b)
        else:
            print("[!] Error: Could not extract valid sequences from the PDB files.")
            
    except Exception as e:
        print(f"\n[!] An error occurred: {e}")
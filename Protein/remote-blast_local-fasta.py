import os
import ssl
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import urllib.error

# Resolve SSL certificate verification issues
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
      ssl._create_default_https_context = _create_unverified_https_context

Entrez.email = "your.email@example.com"

def run_remote_blast(query_sequence, program="blastp", database="nr"):
   
    print(f"\n Running {program.upper()} against '{database}' database on NCBI servers...")
    print(" This may take a few minutes depending on server load...")
    
    try:
        # qblast takes program, database, and query sequence
        result_handle = NCBIWWW.qblast(program, database, query_sequence)
        print(" BLAST search completed successfully.")
        return result_handle
    except urllib.error.URLError as e:
        print(f"\n[!] Network Error during BLAST: {e}")
        return None

def parse_blast_and_save(result_handle, output_filename="significant_hits.fasta", e_value_threshold=0.01):
    print(f"\n Parsing BLAST results (E-value threshold: {e_value_threshold})...")
    blast_records = NCBIXML.parse(result_handle)
    
    saved_count = 0
    with open(output_filename, "w") as out_file:
        # Iterate over all BLAST queries (usually just one if we submitted one sequence)
        for blast_record in blast_records:
            # Iterate over all alignments (hits) for this query
            for alignment in blast_record.alignments:
                # Check the first HSP (High-scoring Segment Pair) in the alignment
                for hsp in alignment.hsps:
                    if hsp.expect < e_value_threshold:
                        # Create FASTA format string
                        header = f">{alignment.title}"
                        # The aligned sequence from the database hit
                        sequence = hsp.sbjct 
                        
                        out_file.write(f"{header}\n{sequence}\n")
                        saved_count += 1
                        
                        
                        print(f" -> Saved Hit: {alignment.title[:50]}... | E-value: {hsp.expect}")
                        break # Move to next alignment after saving the best HSP

    print(f"\n Process complete. Saved {saved_count} significant hits to '{output_filename}'")

if __name__ == "__main__":
    print("--- Remote BLAST & Local FASTA Generator ---")
    
    # A sample protein sequence (shortened for demonstration - this is part of Insulin)
    sample_query = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    
    print("\nDefault query sequence (Human Insulin precursor):")
    print(f"{sample_query[:50]}...")
    
    use_custom = input("\nDo you want to enter your own protein sequence? (y/n): ").strip().lower()
    if use_custom == 'y':
        query_sequence = input("Enter your protein sequence: ").strip()
    else:
        query_sequence = sample_query

    output_file = "blast_results.fasta"
    
    try:
        # Step 1: Run BLAST
        blast_handle = run_remote_blast(query_sequence)
        
        # Step 2: Parse and Save
        if blast_handle:
            parse_blast_and_save(blast_handle, output_filename=output_file)
            blast_handle.close()
            
    except Exception as e:
        print(f"\n[!] An unexpected error occurred: {e}")
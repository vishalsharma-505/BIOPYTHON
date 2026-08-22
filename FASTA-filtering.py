import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analyze_and_filter_fasta(file_path, length_threshold=100):
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    print(f"\nReading and analyzing: {file_path}")
    with open(file_path, "r") as handle:
       
        records = list(SeqIO.parse(handle, "fasta"))
        
    if not records:
        print("[!] File is empty or not in correct FASTA format.")
        return

    lengths = [len(record.seq) for record in records]
    total_seqs = len(records)
    
    # Basic Stats
    print("\n" + "="*30)
    print("      FASTA SUMMARY STATS      ")
    print("="*30)
    print(f"Total Sequences : {total_seqs}")
    print(f"Max Length      : {max(lengths)}")
    print(f"Min Length      : {min(lengths)}")
    print(f"Average Length  : {sum(lengths) / total_seqs:.2f}")
    print("="*30)

    print(f"\nFiltering sequences longer than {length_threshold}...")
    
    filtered_records = []
    
    for record in records:
        seq_length = len(record.seq)
        if seq_length > length_threshold:
            try:
                gc_cont = gc_fraction(record.seq) * 100
                gc_text = f" | GC Content: {gc_cont:.1f}%"
            except Exception:
                gc_text = "" 
                
            print(f" -> KEEP: {record.id} (Length: {seq_length}{gc_text})")
            filtered_records.append(record)
            
    # 5. Exporting (Writing to a new FASTA file)
    if filtered_records:
        output_filename = "filtered_sequences.fasta"
        with open(output_filename, "w") as out_handle:
            SeqIO.write(filtered_records, out_handle, "fasta")
            
        print(f"\nSuccess! Saved {len(filtered_records)} filtered sequences to '{output_filename}'")
    else:
        print(f"\n No sequences found that are longer than {length_threshold}.")

if __name__ == "__main__":
    print("--- Advanced FASTA Data Pipeline ---")
    
    
    target_file = input("Enter the path to your FASTA file (e.g., example.fasta): ").strip()
    
    try:
        min_length = int(input("Enter minimum sequence length to keep: ").strip())
    except ValueError:
        print("Invalid input for length. Using default: 100")
        min_length = 100
        
    analyze_and_filter_fasta(target_file, min_length)
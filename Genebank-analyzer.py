import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analyze_genbank(file_path):
    if not os.path.exists(file_path):
        print(f"Error: Could not find {file_path}")
        return

    print(f"\nProcessing GenBank file: {file_path}")
    
    with open(file_path, "r") as handle:
      
        records = list(SeqIO.parse(handle, "genbank"))
        
    if not records:
        print("No valid GenBank records found.")
        return
        
    for record in records:
        print("\n" + "-"*40)
        print(f"ID: {record.id}")
        print(f"Name: {record.name}")
        print(f"Description: {record.description}")
        print(f"Sequence Length: {len(record.seq)} bp")
        
        try:
            gc = gc_fraction(record.seq) * 100
            print(f"GC Content: {gc:.2f}%")
        except:
            pass
            
        print("\nFeature Summary:")
        # dictionary to count feature types (CDS, gene, etc)
        feature_counts = {}
        coding_sequences = []
        
        for feature in record.features:
            ftype = feature.type
            feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
            
            if ftype == "CDS":
               try:
                   product = feature.qualifiers.get('product', ['Unknown'])[0]
                   protein_id = feature.qualifiers.get('protein_id', ['None'])[0]
                   coding_sequences.append(f"{product} ({protein_id})")
               except Exception:
                   continue

        for f_type, count in feature_counts.items():
            print(f" - {f_type}: {count}")
        if coding_sequences:
            print(f"\nExtracted {len(coding_sequences)} Coding Sequences (CDS):")
            for i, cds in enumerate(coding_sequences[:5]): # show first 5
                print(f" {i+1}. {cds}")
            if len(coding_sequences) > 5:
                print(f" ... and {len(coding_sequences)-5} more.")
            

if __name__ == "__main__":
   filepath = input("Enter path to GenBank (.gb or .gbk) file: ").strip()
   analyze_genbank(filepath)
   
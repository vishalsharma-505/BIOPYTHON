import os
from Bio.Blast import NCBIXML

def analyze_blast_xml(xml_file):
    if not os.path.exists(xml_file):
        print("file not found check path")
        return

    with open(xml_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        for record in blast_records:
            print(f"Query: {record.query}")
            
            total_hits = len(record.alignments)
            if total_hits == 0:
                print("no hits found in this xml")
                continue
                
            identities = []
            scores = []
            evalues = []
            
            print("\n--- Top Pairwise Alignments ---")
            
            for i, alignment in enumerate(record.alignments):
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100
                    identities.append(identity_pct)
                    scores.append(hsp.score)
                    evalues.append(hsp.expect)
                    
                    if i < 3:
                        print(f"\nHit: {alignment.title[:60]}...")
                        print(f"E-value: {hsp.expect} | Identity: {identity_pct:.2f}%")
                        print("Query :", hsp.query[0:60])
                        print("Match :", hsp.match[0:60])
                        print("Sbjct :", hsp.sbjct[0:60])
            
            print("\n--- Summary Stats ---")
            print(f"Total Hits Analysed : {total_hits}")
            print(f"Highest Score       : {max(scores)}")
            print(f"Best E-value        : {min(evalues)}")
            print(f"Average Identity    : {sum(identities)/len(identities):.2f}%")

if __name__ == "__main__":
    filepath = input("Enter BLAST XML file path: ")
    analyze_blast_xml(filepath)
from Bio import Entrez
import ssl
import urllib.error

# Resolves SSL certificate verification issues on some systems
ssl._create_default_https_context = ssl._create_unverified_context

Entrez.email = "your.email@example.com" 

def fetch_entrez_data(database, search_term, batch_size=50):
    """
    Searches the NCBI database and allows the user to interactively 
    fetch FASTA sequences in batches.
    """
    start = 0
    
    print(f"\nSearching '{database}' database for '{search_term}'...")
    
    while True:
        try:
            # 1. Search NCBI
            handle = Entrez.esearch(db=database, term=search_term, retmax=batch_size, retstart=start)
            results = Entrez.read(handle)
            handle.close() 
            
            total_records = int(results["Count"])
            id_list = results["IdList"]
            
            if not id_list:
                print("No records found or end of results reached.")
                break
                
            print(f"\n--- Showing results {start + 1} to {start + len(id_list)} out of {total_records} ---")
            for ncbi_id in id_list:
                print(ncbi_id)
                
            # 2. User Interaction
            print("\nOptions:")
            print("- Enter an [ID Number] to fetch its FASTA sequence")
            print("- Enter [n] to load the next batch")
            print("- Enter [q] to quit")
            
            user_choice = input("Your choice: ").strip().lower()
            
            if user_choice == 'q':
                print("Exiting program.")
                break
                
            elif user_choice == 'n':
                start += batch_size
                if start >= total_records:
                    print("\nAll records have been viewed.")
                    break
                    
            elif user_choice in id_list:
                print(f"\nFetching FASTA data for ID: {user_choice}...\n")


                # 3. Fetch Data
                fetch_handle = Entrez.efetch(db=database, id=user_choice, rettype="fasta", retmode="text")
                print(fetch_handle.read())
                fetch_handle.close()
                
                # Move to next batch after viewing a sequence
                start += batch_size 
            else:
                print("\nInvalid input. Please enter a valid ID from the list, 'n', or 'q'.")
                
        except urllib.error.URLError as e:
            print(f"\nNetwork Error: {e}. Please check your internet connection.")
            break
        except Exception as e:
            print(f"\nAn unexpected error occurred: {e}")
            break

if __name__ == "__main__":
    print("=== NCBI Entrez Sequence Fetcher ===")
    user_db = input("Enter database (e.g., nucleotide, protein): ").strip().lower()
    user_term = input("Enter search term (e.g., insulin, cancer): ").strip()
    
    fetch_entrez_data(user_db, user_term)
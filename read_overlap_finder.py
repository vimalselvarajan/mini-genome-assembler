from Bio import SeqIO
from collections import defaultdict

def read_fasta(file_path):
    reads = {}
    for record in SeqIO.parse(file_path, "fasta"):
        reads[record.id] = str(record.seq)
    return reads

def build_kmer_index(reads, k=40):
    kmer_index = defaultdict(set)
    for read_id, seq in reads.items():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_index[kmer].add(read_id)
    return kmer_index

def find_best_buddies(reads, k=40, min_overlap=40):
    kmer_index = build_kmer_index(reads, k)
    best_buddies = {}

    for A_id, A_seq in reads.items():
        max_overlap = 0
        best_B = None
        tie = False

        suffix_kmer = A_seq[-k:]

        candidates = kmer_index.get(suffix_kmer, set()) - {A_id}

        for B_id in candidates:
            B_seq = reads[B_id]
            max_possible = min(len(A_seq), len(B_seq))
            olen = 0

            for length in range(max_possible, min_overlap - 1, -1):
                if A_seq[-length:] == B_seq[:length]:
                    olen = length
                    break

            if olen > max_overlap:
                max_overlap = olen
                best_B = B_id
                tie = False
                
            elif olen == max_overlap and olen >= min_overlap:
                tie = True

        if max_overlap >= min_overlap and not tie:
            best_buddies[A_id] = (best_B, max_overlap)

    return best_buddies

if __name__ == "__main__":
    fasta_file = "sample_data/reads.fa"  
    reads = read_fasta(fasta_file)
    buddies = find_best_buddies(reads)
    
    with open("overlaps.txt", "w") as f:
        for A, (B, olen) in buddies.items():
            f.write(f"{A} {B} {olen}\n")

    print(f"Done. Wrote {len(buddies)} overlaps to overlaps.txt")
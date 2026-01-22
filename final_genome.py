from Bio import SeqIO

def load_reads(fasta_file):
    read_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        read_dict[record.id] = str(record.seq)
    return read_dict

def parse_unitigs_with_overlaps(unitig_file):
    unitigs = []
    with open(unitig_file) as f:
        current_ids = []
        current_ovls = []
        in_unitig = False
        for line in f:
            if line.startswith('START OF UNITIG'):
                in_unitig = True
                parts = line.strip().split()
                first_id = parts[-1]
                current_ids = [first_id]
                current_ovls = []
                continue
            elif line.startswith('END OF UNITIG'):
                in_unitig = False
                if current_ids:
                    unitigs.append((current_ids, current_ovls))
                continue
            if in_unitig:
                parts = line.strip().split()
                if len(parts) == 2:
                    current_ids.append(parts[0])
                    current_ovls.append(int(parts[1]))
    return unitigs

def build_unitig_sequence(read_ids, ovls, read_dict):
    seq = read_dict[read_ids[0]]
    for i in range(1, len(read_ids)):
        ovl = ovls[i-1]  # ovls aligns with transition from i-1 to i
        seq += read_dict[read_ids[i]][ovl:]
    return seq

def write_unitigs_to_fasta(unitigs, read_dict, filename):
    with open(filename, 'w') as f:
        for i, (read_ids, ovls) in enumerate(unitigs):
            seq = build_unitig_sequence(read_ids, ovls, read_dict)
            f.write(f'>unitig_{i}\n')
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')

# ---- Main workflow ----

reads_fasta = "sample_data/reads.fa"
unitigs_txt = "unitigs.txt"
output_unitigs_fasta = "unitigs.fasta"

read_dict = load_reads(reads_fasta)
unitigs = parse_unitigs_with_overlaps(unitigs_txt)
write_unitigs_to_fasta(unitigs, read_dict, output_unitigs_fasta)

print(f"Done! Wrote {len(unitigs)} unitigs to {output_unitigs_fasta}")

# Mini Genome Assembler

## Project Description

An educational mini genome assembler for **CS 144 (Algorithms for Bioinformatics)** that identifies maximum suffix–prefix overlaps, builds unitigs, and reconstructs a genome from error-free reads.

This project implements a simplified genome assembly pipeline using **exact string overlaps**. Given a set of error-free sequencing reads from the forward strand, the program:

1. Identifies **maximum suffix–prefix overlaps** between reads using a **k-mer–based indexing strategy**
2. Constructs **unitigs** by joining reads that are **mutual maximum-overlap neighbors**
3. Assembles a **final genome sequence** by overlapping unitigs

The goal of this project is to illustrate core algorithmic ideas behind **Overlap–Layout–Consensus (OLC)** assembly.

---

## Requirements

- Python 3.8+
- BioPython
- NetworkX

Install dependencies:
```bash
pip install biopython networkx
```
---

## Commands

Find maximum suffix–prefix overlaps
```bash
python read_overlap_finder.py
```
Build unitigs
```bash
python assemble_unitigs.py
```
Assemble the final genome
```bash
python final_genome.py
```

## Pipeline Structure

### Step 1: Maximum Suffix–Prefix Overlaps
- For each read **A**, identify the read **B** with the **longest exact suffix–prefix overlap**
- Overlaps must be **at least 40 nucleotides**
- If the longest overlap is not unique, the read has no valid outgoing overlap
- Implemented efficiently using **k-mer indexing** (k = 40)

**Output:** `overlaps.txt`  

---

### Step 2: Unitig Construction
- Builds a directed overlap graph from the overlap file
- Joins reads into **unitigs** when they are **mutual maximum-overlap neighbors**
- Each read appears in exactly one unitig
- Implemented using **NetworkX**

**Output:** `unitigs.txt`

---

### Step 3: Final Genome Assembly
- Converts unitigs into sequences
- Determines the correct ordering by overlapping unitigs
- Produces the final assembled genome in FASTA format
- Unitigs may be reused if necessary

**Output:** `unitigs.fasta` / final genome FASTA  
(Expected genome length: **7959 bp**)





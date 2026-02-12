# app.py

import streamlit as st
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import numpy as np

st.set_page_config(page_title="Forensic DNA Platform", layout="wide")

st.title("🧬 Advanced Forensic DNA Analysis Platform")

# ================= FUNCTIONS =================


def dna_calculator(dna):
    dna = dna.upper()
    a = t = g = c = 0

    for base in dna:
        if base == "A":
            a += 1
        elif base == "T":
            t += 1
        elif base == "G":
            g += 1
        elif base == "C":
            c += 1

    length = a + t + g + c
    gc_content = (g + c) / length * 100 if length > 0 else 0
    return length, a, t, g, c, gc_content


def find_motif(dna, motif):
    dna = dna.upper()
    motif = motif.upper()
    positions = []

    for i in range(len(dna) - len(motif) + 1):
        if dna[i:i+len(motif)] == motif:
            positions.append(i + 1)

    return positions


def count_str_repeats(sequence, motif):
    sequence = sequence.upper()
    motif = motif.upper()

    if motif == "":
        return 0, 0

    repeat_count = sequence.count(motif)
    allele_size = repeat_count * len(motif)

    return repeat_count, allele_size


def compare_sequences(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    differences = []
    min_length = min(len(seq1), len(seq2))

    for i in range(min_length):
        if seq1[i] != seq2[i]:
            differences.append((i+1, seq1[i], seq2[i]))

    return differences


def sequence_similarity(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    min_len = min(len(seq1), len(seq2))
    matches = 0

    for i in range(min_len):
        if seq1[i] == seq2[i]:
            matches += 1

    similarity = (matches / min_len) * 100 if min_len > 0 else 0
    return similarity


def sequence_distance(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    min_len = min(len(seq1), len(seq2))
    diff = 0

    for i in range(min_len):
        if seq1[i] != seq2[i]:
            diff += 1

    return diff


# ================= DNA ANALYSIS =================

st.header("🔬 DNA Sequence Analysis")

sequence_input = st.text_area(
    "Enter DNA sequences (separate multiple with commas)")
motif = st.text_input("Motif search (optional)")

if st.button("Analyze DNA"):
    sequences = sequence_input.split(",")

    for idx, seq in enumerate(sequences, 1):
        seq = seq.strip()
        length, a, t, g, c, gc = dna_calculator(seq)

        st.subheader(f"Sequence {idx}")
        st.write("Length:", length)
        st.write(f"A: {a} | T: {t} | G: {g} | C: {c}")
        st.write(f"GC Content: {gc:.2f}%")

        fig, ax = plt.subplots()
        ax.bar(["A", "T", "G", "C"], [a, t, g, c])
        ax.set_title("Base Count")
        st.pyplot(fig)

        if motif:
            positions = find_motif(seq, motif)
            if positions:
                st.write("Motif found at positions:", positions)
            else:
                st.write("Motif not found.")

st.divider()

# ================= STR =================

st.header("🧬 STR Repeat Analysis")

str_sequence = st.text_input("Enter STR sequence")
str_motif = st.text_input("Enter STR motif (e.g., AGAT)")

if st.button("Analyze STR"):
    repeats, allele = count_str_repeats(str_sequence, str_motif)

    st.write("Repeat Count:", repeats)
    st.write("Estimated Allele Size (bp):", allele)

    fig, ax = plt.subplots()
    ax.bar([str_motif.upper()], [repeats])
    ax.set_title("STR Repeat Count")
    st.pyplot(fig)

st.divider()

# ================= SNP =================

st.header("🧬 SNP Comparator")

seq1 = st.text_input("Enter Sequence 1")
seq2 = st.text_input("Enter Sequence 2")

if st.button("Compare Sequences"):
    snps = compare_sequences(seq1, seq2)

    if snps:
        st.write("SNP Differences Found:")
        for pos, b1, b2 in snps:
            st.write(f"Position {pos}: {b1} → {b2}")

        fig, ax = plt.subplots()
        ax.bar(["SNP Count"], [len(snps)])
        ax.set_title("Number of SNP Differences")
        st.pyplot(fig)
    else:
        st.write("No SNP differences found.")

st.divider()

# ================= BOLD-LIKE SYSTEM =================

st.header("🌍 DNA Barcoding Identification")

reference_db = {
    "Human": "ATGCGTACGTAGCTAGCTAG",
    "Dog": "ATGCGTTCGTAGCTAGTTAG",
    "Cat": "ATGCGTACGTAGATAGCTAG"
}

unknown_seq = st.text_input("Enter Unknown Sample Sequence")

if st.button("Identify Species"):
    results = {}

    for species, ref_seq in reference_db.items():
        similarity = sequence_similarity(unknown_seq, ref_seq)
        results[species] = similarity

    best_match = max(results, key=results.get)

    st.write("Similarity Results:")
    for species, sim in results.items():
        st.write(f"{species}: {sim:.2f}%")

    st.success(f"Most Likely Match: {best_match}")

# ================= PHYLOGENETIC TREE =================

st.header("🌳 Phylogenetic Tree Builder")

tree_input = st.text_area("Enter sequences for tree (separate with commas)")

if st.button("Build Phylogenetic Tree"):

    sequences = [seq.strip() for seq in tree_input.split(",") if seq.strip()]

    if len(sequences) < 2:
        st.warning("Enter at least 2 sequences.")
    else:
        labels = [f"Seq{i+1}" for i in range(len(sequences))]

        # Build condensed distance matrix
        distances = []
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                distances.append(sequence_distance(sequences[i], sequences[j]))

        linked = linkage(distances, method='single')

        fig, ax = plt.subplots()
        dendrogram(linked, labels=labels)
        plt.title("Phylogenetic Tree")
        st.pyplot(fig)

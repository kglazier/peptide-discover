"""File format read/write helpers."""

from pathlib import Path


def read_fasta(path: Path) -> str:
    """Read first sequence from a FASTA file."""
    from Bio import SeqIO

    record = next(SeqIO.parse(str(path), "fasta"))
    return str(record.seq)


def read_pdb_sequence(path: Path) -> str:
    """Extract amino acid sequence from a PDB file."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import protein_letters_3to1

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", str(path))
    model = structure[0]

    sequences = []
    for chain in model.get_chains():
        seq = []
        for residue in chain.get_residues():
            resname = residue.get_resname()
            if resname in protein_letters_3to1:
                seq.append(protein_letters_3to1[resname])
        if seq:
            sequences.append("".join(seq))

    if not sequences:
        raise ValueError(f"No protein chains found in {path}")
    return sequences[0]  # Return first chain


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write sequences to a FASTA file."""
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")

#!/usr/bin/env python3
"""Convert FASTA format files to AlphaFold 3 input JSON format.

This script converts protein/RNA/DNA FASTA files to the AlphaFold 3 input JSON format,
supporting single or multiple sequences with automatic chain assignment.

Usage:
    python fasta_to_af3.py input.fasta -o output.json
    python fasta_to_af3.py input.fasta --mol-type protein -o output.json
    python fasta_to_af3.py input.fasta --mol-type rna -o output.json
    python fasta_to_af3.py input.fasta --name my_protein --seeds 1,2,3 -o output.json
"""

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


class MoleculeType(Enum):
    """Supported molecule types for AlphaFold 3 input."""
    PROTEIN = "protein"
    RNA = "rna"
    DNA = "dna"


@dataclass
class FASTAEntry:
    """Represents a single FASTA entry."""
    header: str
    sequence: str
    identifier: str = ""
    description: str = ""

    def __post_init__(self):
        """Parse header into identifier and description."""
        if not self.identifier:
            parts = self.header.lstrip(">").split(None, 1)
            self.identifier = parts[0] if parts else ""
            self.description = parts[1] if len(parts) > 1 else ""


def parse_fasta(fasta_path: str) -> list[FASTAEntry]:
    """Parse a FASTA file and return list of FASTAEntry objects.
    
    Args:
        fasta_path: Path to FASTA file.
        
    Returns:
        List of FASTAEntry objects.
        
    Raises:
        FileNotFoundError: If FASTA file doesn't exist.
        ValueError: If FASTA file is empty or invalid.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    with open(fasta_path, 'r') as f:
        content = f.read().strip()
    
    if not content:
        raise ValueError(f"FASTA file is empty: {fasta_path}")
    
    entries = []
    current_header = None
    current_sequence = []
    
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('>'):
            # Save previous entry if exists
            if current_header is not None:
                sequence = ''.join(current_sequence).upper()
                if sequence:
                    entries.append(FASTAEntry(
                        header=current_header,
                        sequence=sequence,
                    ))
            
            current_header = line[1:]  # Remove '>'
            current_sequence = []
        else:
            # Remove any whitespace and numbers from sequence
            clean_line = re.sub(r'[\s\d]', '', line)
            current_sequence.append(clean_line)
    
    # Save last entry
    if current_header is not None:
        sequence = ''.join(current_sequence).upper()
        if sequence:
            entries.append(FASTAEntry(
                header=current_header,
                sequence=sequence,
            ))
    
    if not entries:
        raise ValueError(f"No valid FAST entries found in: {fasta_path}")
    
    return entries


def validate_protein_sequence(sequence: str) -> tuple[bool, list[str]]:
    """Validate a protein sequence.
    
    Args:
        sequence: Protein sequence string.
        
    Returns:
        Tuple of (is_valid, list_of_warnings).
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    warnings = []
    
    invalid_chars = set(sequence) - valid_aa - {'X', 'B', 'Z', 'U', 'O', '*'}
    if invalid_chars:
        return False, [f"Invalid amino acid characters found: {invalid_chars}"]
    
    if 'X' in sequence:
        warnings.append("Sequence contains unknown amino acid 'X'")
    
    if '*' in sequence:
        warnings.append("Sequence contains stop codon '*'")
    
    return True, warnings


def validate_nucleic_acid_sequence(sequence: str, mol_type: str) -> tuple[bool, list[str]]:
    """Validate a nucleic acid sequence.
    
    Args:
        sequence: Nucleic acid sequence string.
        mol_type: Either 'rna' or 'dna'.
        
    Returns:
        Tuple of (is_valid, list_of_warnings).
    """
    if mol_type == "rna":
        valid_bases = set('ACGU')
    else:  # dna
        valid_bases = set('ACGT')
    
    warnings = []
    invalid_chars = set(sequence) - valid_bases - {'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}
    
    if invalid_chars:
        return False, [f"Invalid nucleotide characters found: {invalid_chars}"]
    
    return True, warnings


def assign_chains(num_entries: int) -> list[str]:
    """Assign chain IDs for multiple sequences.
    
    Args:
        num_entries: Number of sequences.
        
    Returns:
        List of chain IDs (A, B, C, ..., Z, AA, AB, ...).
    """
    chains = []
    for i in range(num_entries):
        # Single letter: A-Z
        if i < 26:
            chains.append(chr(ord('A') + i))
        # Double letters: AA, AB, ...
        else:
            first = chr(ord('A') + (i // 26) - 1)
            second = chr(ord('A') + (i % 26))
            chains.append(first + second)
    return chains


def detect_sequence_type(sequence: str, identifier: str = "", verbose: bool = False) -> MoleculeType:
    """Detect the molecule type of a single sequence.
    
    Analyzes individual sequence to determine if it's protein, RNA, or DNA.
    
    Args:
        sequence: The sequence string.
        identifier: Sequence identifier (for verbose output).
        verbose: Whether to print detection info.
        
    Returns:
        Detected MoleculeType for this sequence.
    """
    seq = sequence.upper()
    seq_len = len(seq)
    
    if seq_len == 0:
        return MoleculeType.PROTEIN
    
    # Count key characters
    has_u = 'U' in seq
    has_t = 'T' in seq
    
    # Protein-specific residues (not in nucleic acids)
    protein_specific = set('DEFILMPQVWY')  # Definitive amino acids
    protein_count = sum(1 for c in seq if c in protein_specific)
    
    # Nucleic acid characters
    nucleic_chars = set('ACGTU')
    nucleic_count = sum(1 for c in seq if c in nucleic_chars)
    nucleic_ratio = nucleic_count / seq_len
    
    # Detection logic
    # Rule 1: If has protein-specific residues, it's protein
    if protein_count > 0:
        if verbose:
            print(f"  {identifier}: Protein ({protein_count} protein-specific residues)")
        return MoleculeType.PROTEIN
    
    # Rule 2: If mostly nucleic acids, check U vs T
    if nucleic_ratio > 0.8:
        if has_u and not has_t:
            if verbose:
                print(f"  {identifier}: RNA (has U, no T)")
            return MoleculeType.RNA
        elif has_t and not has_u:
            if verbose:
                print(f"  {identifier}: DNA (has T, no U)")
            return MoleculeType.DNA
        elif has_u and has_t:
            # Both present - likely protein with unusual composition
            if verbose:
                print(f"  {identifier}: Protein (both U and T present)")
            return MoleculeType.PROTEIN
        else:
            # Only A, C, G - ambiguous, default to protein
            if verbose:
                print(f"  {identifier}: Protein (only ACG, defaulting)")
            return MoleculeType.PROTEIN
    
    # Rule 3: Mixed composition - likely protein
    if verbose:
        print(f"  {identifier}: Protein (mixed composition)")
    return MoleculeType.PROTEIN


def create_af3_input_mixed(
    fasta_entries: list[FASTAEntry],
    name: str,
    seeds: list[int],
    verbose: bool = False,
    version: int = 1,
) -> dict:
    """Create AlphaFold 3 input JSON with mixed molecule types.
    
    Each sequence is analyzed individually and assigned the appropriate
    molecule type (protein, rna, or dna).
    
    Args:
        fasta_entries: List of FASTAEntry objects.
        name: Target name for the complex.
        seeds: List of random seeds.
        verbose: Whether to print detection info.
        version: AlphaFold 3 input version.
        
    Returns:
        Dictionary representing AlphaFold 3 input JSON.
    """
    chains = assign_chains(len(fasta_entries))
    
    sequences = []
    type_counts = {MoleculeType.PROTEIN: 0, MoleculeType.RNA: 0, MoleculeType.DNA: 0}
    
    for entry, chain in zip(fasta_entries, chains):
        # Detect type for this specific sequence
        mol_type = detect_sequence_type(entry.sequence, entry.identifier, verbose)
        type_counts[mol_type] += 1
        
        if verbose:
            print(f"  Chain {chain}: {mol_type.value} ({len(entry.sequence)} residues)")
        
        # Create sequence entry
        seq_entry = {
            mol_type.value: {
                "id": [chain],
                "sequence": entry.sequence,
            }
        }
        sequences.append(seq_entry)
    
    if verbose:
        print(f"\nMolecule type summary:")
        for mol_type, count in type_counts.items():
            if count > 0:
                print(f"  {mol_type.value}: {count} chain(s)")
    
    af3_input = {
        "name": name,
        "sequences": sequences,
        "modelSeeds": seeds,
        "dialect": "alphafold3",
        "version": version,
    }
    
    return af3_input


def create_af3_input(
    fasta_entries: list[FASTAEntry],
    name: str,
    mol_type: MoleculeType,
    seeds: list[int],
    version: int = 1,
) -> dict:
    """Create AlphaFold 3 input JSON structure.
    
    Args:
        fasta_entries: List of FASTAEntry objects.
        name: Target name for the complex.
        mol_type: Molecule type (protein, rna, or dna).
        seeds: List of random seeds.
        version: AlphaFold 3 input version.
        
    Returns:
        Dictionary representing AlphaFold 3 input JSON.
    """
    chains = assign_chains(len(fasta_entries))
    
    sequences = []
    for entry, chain in zip(fasta_entries, chains):
        if mol_type == MoleculeType.PROTEIN:
            seq_entry = {
                mol_type.value: {
                    "id": [chain],
                    "sequence": entry.sequence,
                }
            }
        elif mol_type in (MoleculeType.RNA, MoleculeType.DNA):
            seq_entry = {
                mol_type.value: {
                    "id": [chain],
                    "sequence": entry.sequence,
                }
            }
        else:
            raise ValueError(f"Unsupported molecule type: {mol_type}")
        
        sequences.append(seq_entry)
    
    af3_input = {
        "name": name,
        "sequences": sequences,
        "modelSeeds": seeds,
        "dialect": "alphafold3",
        "version": version,
    }
    
    return af3_input


def generate_name(fasta_entries: list[FASTAEntry], fasta_path: str) -> str:
    """Generate a target name from FASTA entries or file path.
    
    Args:
        fasta_entries: List of FASTAEntry objects.
        fasta_path: Path to FASTA file.
        
    Returns:
        Generated target name.
    """
    # Try to use first entry's identifier
    if fasta_entries and fasta_entries[0].identifier:
        name = fasta_entries[0].identifier
        # Clean up identifier
        name = re.sub(r'[^a-zA-Z0-9_-]', '_', name)
        return name[:50]  # Truncate if too long
    
    # Fall back to filename
    return Path(fasta_path).stem


def detect_molecule_type(
    fasta_entries: list[FASTAEntry],
    verbose: bool = False,
) -> tuple[MoleculeType, list[str]]:
    """Detect molecule type from sequence content using comprehensive analysis.
    
    Uses multiple heuristics to determine if sequences are protein, RNA, or DNA:
    1. Character composition analysis
    2. Sequence length patterns
    3. Amino acid vs nucleotide frequency
    4. Scoring-based classification
    
    Args:
        fasta_entries: List of FASTAEntry objects.
        verbose: Whether to print detailed detection info.
        
    Returns:
        Tuple of (detected MoleculeType, list of explanation strings).
    """
    explanations = []
    
    # Collect all sequences
    all_sequences = [entry.sequence.upper() for entry in fasta_entries]
    full_sequence = ''.join(all_sequences)
    seq_length = len(full_sequence)
    
    if seq_length == 0:
        return MoleculeType.PROTEIN, ["Empty sequence, defaulting to protein"]
    
    # Count character frequencies
    counts = {
        'A': full_sequence.count('A'),
        'C': full_sequence.count('C'),
        'G': full_sequence.count('G'),
        'T': full_sequence.count('T'),
        'U': full_sequence.count('U'),
    }
    
    # Amino acids that are NOT in nucleic acid bases
    # Nucleic acids only use: A, C, G, T, U (and ambiguity codes)
    # These amino acids definitively indicate protein
    protein_specific = {
        'D': full_sequence.count('D'),  # Aspartic acid
        'E': full_sequence.count('E'),  # Glutamic acid
        'F': full_sequence.count('F'),  # Phenylalanine
        'I': full_sequence.count('I'),  # Isoleucine
        'L': full_sequence.count('L'),  # Leucine
        'M': full_sequence.count('M'),  # Methionine
        'N': full_sequence.count('N'),  # Asparagine
        'P': full_sequence.count('P'),  # Proline
        'Q': full_sequence.count('Q'),  # Glutamine
        'R': full_sequence.count('R'),  # Arginine
        'S': full_sequence.count('S'),  # Serine
        'V': full_sequence.count('V'),  # Valine
        'W': full_sequence.count('W'),  # Tryptophan
        'Y': full_sequence.count('Y'),  # Tyrosine
        'K': full_sequence.count('K'),  # Lysine
        'H': full_sequence.count('H'),  # Histidine
    }
    
    # Ambiguity codes (can be in both protein and nucleic acid)
    # Protein: B (Asx), Z (Glx), X (unknown)
    # Nucleic: N, R, Y, S, W, K, M, B, D, H, V
    protein_ambiguous = {
        'B': full_sequence.count('B'),  # Asx (Asp or Asn) or nucleic acid
        'Z': full_sequence.count('Z'),  # Glx (Glu or Gln)
        'X': full_sequence.count('X'),  # Unknown amino acid
        'J': full_sequence.count('J'),  # Leu or Ile
        'O': full_sequence.count('O'),  # Pyrrolysine
        '*': full_sequence.count('*'),  # Stop codon
    }
    
    nucleic_ambiguous = {
        'N': full_sequence.count('N'),  # Any nucleotide
        'R': full_sequence.count('R'),  # Purine (A or G)
        'Y': full_sequence.count('Y'),  # Pyrimidine (C or T/U)
        'S': full_sequence.count('S'),  # Strong (C or G)
        'W': full_sequence.count('W'),  # Weak (A or T/U)
        'K': full_sequence.count('K'),  # Keto (G or T/U)
        'M': full_sequence.count('M'),  # Amino (A or C)
        'D': full_sequence.count('D'),  # Not C (A, G, or T/U)
        'H': full_sequence.count('H'),  # Not G (A, C, or T/U)
        'V': full_sequence.count('V'),  # Not T/U (A, C, or G)
    }
    
    # Calculate scores
    protein_score = 0
    nucleic_score = 0
    
    # Rule 1: Presence of protein-specific residues (strongest indicator)
    total_protein_specific = sum(protein_specific.values())
    if total_protein_specific > 0:
        protein_score += total_protein_specific * 10
        specific_found = [aa for aa, count in protein_specific.items() if count > 0]
        explanations.append(
            f"Protein-specific residues found: {specific_found} "
            f"({total_protein_specific} total, strong protein indicator)"
        )
    
    # Rule 2: Check for U vs T (nucleic acid discrimination)
    if counts['U'] > 0 and counts['T'] == 0:
        nucleic_score += counts['U'] * 5 + 50  # Bonus for U without T
        explanations.append(
            f"Uracil (U) present without Thymine (T): strong RNA indicator"
        )
    elif counts['T'] > 0 and counts['U'] == 0:
        nucleic_score += counts['T'] * 5 + 30  # Bonus for T without U
        explanations.append(
            f"Thymine (T) present without Uracil (U): DNA indicator"
        )
    elif counts['U'] > 0 and counts['T'] > 0:
        # Both present - likely protein with ambiguous codes
        explanations.append(
            f"Both U ({counts['U']}) and T ({counts['T']}) present - likely protein"
        )
        protein_score += 20
    
    # Rule 3: Amino acid composition analysis
    # Proteins typically have diverse amino acid distribution
    # Nucleic acids are dominated by A, C, G, T/U
    valid_protein = set('ACDEFGHIKLMNPQRSTVWY')
    valid_nucleic = set('ACGTUNRYSWKMBDHV')
    
    protein_chars = sum(1 for c in full_sequence if c in valid_protein)
    nucleic_chars = sum(1 for c in full_sequence if c in valid_nucleic)
    
    protein_ratio = protein_chars / seq_length
    nucleic_ratio = nucleic_chars / seq_length
    
    if protein_ratio > 0.8 and protein_ratio > nucleic_ratio:
        protein_score += int(protein_ratio * 30)
        explanations.append(
            f"High protein character ratio: {protein_ratio:.1%} "
            f"(supports protein classification)"
        )
    elif nucleic_ratio > 0.9 and nucleic_ratio > protein_ratio:
        nucleic_score += int(nucleic_ratio * 20)
        explanations.append(
            f"High nucleic acid character ratio: {nucleic_ratio:.1%} "
            f"(supports nucleic acid classification)"
        )
    
    # Rule 4: Sequence length patterns
    # Very long sequences (>1000) are more likely proteins
    # Short sequences (<100) are more likely nucleic acids
    avg_length = sum(len(seq) for seq in all_sequences) / len(all_sequences)
    if avg_length > 500:
        protein_score += 15
        explanations.append(
            f"Average sequence length {avg_length:.0f} residues "
            f"(typical for proteins)"
        )
    elif avg_length < 100:
        nucleic_score += 10
        explanations.append(
            f"Average sequence length {avg_length:.0f} residues "
            f"(could be nucleic acid)"
        )
    
    # Rule 5: Check for amino acids that are rare in nucleic acid ambiguity
    # D is particularly telling because it's not in nucleic codes
    if protein_specific['D'] > 0 or protein_specific['E'] > 0:
        protein_score += 25
        explanations.append(
            f"Acidic residues D/E present (definitive protein markers)"
        )
    
    # Rule 6: Analyze if ambiguous codes are protein or nucleic
    # In proteins, B/Z/X are common; in nucleic acids, N/R/Y/S are common
    total_prot_amBIG = sum(protein_ambiguous.values())
    total_nuc_amBIG = sum(nucleic_ambiguous.values())
    
    if total_prot_amBIG > total_nuc_amBIG and total_prot_amBIG > 0:
        protein_score += total_prot_amBIG * 3
        explanations.append(
            f"Protein ambiguity codes (B/Z/X) more prevalent than "
            f"nucleic ambiguity codes"
        )
    elif total_nuc_amBIG > total_prot_amBIG and total_nuc_amBIG > 0:
        nucleic_score += total_nuc_amBIG * 2
        explanations.append(
            f"Nucleic acid ambiguity codes (N/R/Y/S) more prevalent"
        )
    
    # Final decision
    explanations.append(f"\nFinal scores: Protein={protein_score}, Nucleic={nucleic_score}")
    
    if protein_score > nucleic_score:
        if verbose:
            print("Auto-detection analysis:")
            for exp in explanations:
                print(f"  {exp}")
        return MoleculeType.PROTEIN, explanations
    
    # Distinguish RNA vs DNA
    if counts['U'] > counts['T']:
        if verbose:
            print("Auto-detection analysis:")
            for exp in explanations:
                print(f"  {exp}")
        return MoleculeType.RNA, explanations
    else:
        if verbose:
            print("Auto-detection analysis:")
            for exp in explanations:
                print(f"  {exp}")
        return MoleculeType.DNA, explanations


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Convert FASTA format files to AlphaFold 3 input JSON format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (auto-detect molecule type)
  python fasta_to_af3.py input.fasta -o output.json

  # Specify molecule type
  python fasta_to_af3.py input.fasta --mol-type protein -o output.json
  python fasta_to_af3.py input.fasta --mol-type rna -o output.json
  python fasta_to_af3.py input.fasta --mol-type dna -o output.json

  # Custom name and seeds
  python fasta_to_af3.py input.fasta --name my_complex --seeds 1,2,3 -o output.json

  # Multiple sequences (will be assigned chains A, B, C, ...)
  python fasta_to_af3.py multi.fasta --mol-type protein -o output.json

  # Verbose mode with validation
  python fasta_to_af3.py input.fasta -o output.json --verbose
        """,
    )
    
    parser.add_argument(
        'fasta',
        help='Input FASTA file path',
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output JSON file path',
    )
    parser.add_argument(
        '--mol-type',
        choices=['protein', 'rna', 'dna'],
        default=None,
        help='Molecule type (default: auto-detect from sequence content)',
    )
    parser.add_argument(
        '--name',
        default=None,
        help='Target name (default: auto-generate from first entry or filename)',
    )
    parser.add_argument(
        '--seeds',
        default='1',
        help='Comma-separated list of random seeds (default: 1)',
    )
    parser.add_argument(
        '--version',
        type=int,
        default=1,
        help='AlphaFold 3 input version (default: 1)',
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output',
    )
    parser.add_argument(
        '--no-validate',
        action='store_true',
        help='Disable sequence validation',
    )
    
    args = parser.parse_args()
    
    # Parse seeds
    try:
        seeds = [int(s.strip()) for s in args.seeds.split(',')]
        if any(s < 0 for s in seeds):
            raise ValueError("Seeds must be non-negative integers")
    except ValueError as e:
        print(f"Error: Invalid seeds argument '{args.seeds}': {e}", file=sys.stderr)
        sys.exit(1)
    
    # Parse FASTA file
    if args.verbose:
        print(f"Parsing FASTA file: {args.fasta}")
    
    try:
        fasta_entries = parse_fasta(args.fasta)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print(f"Found {len(fasta_entries)} sequence(s)")
        for i, entry in enumerate(fasta_entries):
            print(f"  {i+1}. {entry.identifier} ({len(entry.sequence)} residues)")

    # Auto-detect molecule type if not specified
    if args.mol_type is None:
        if args.verbose:
            print("Auto-detect: Using mixed mode (each sequence analyzed individually)")
        
        # Use mixed mode - detect type for each sequence individually
        mol_type = None  # Will be handled by create_af3_input_mixed
        use_mixed_mode = True
    else:
        mol_type = MoleculeType(args.mol_type)
        use_mixed_mode = False
    
    if args.verbose and mol_type:
        print(f"Using molecule type: {mol_type.value}")
    
    # Validate sequences
    if not args.no_validate and not use_mixed_mode:
        if args.verbose:
            print("Validating sequences...")
        
        all_valid = True
        for i, entry in enumerate(fasta_entries):
            if mol_type == MoleculeType.PROTEIN:
                is_valid, warnings = validate_protein_sequence(entry.sequence)
            else:
                is_valid, warnings = validate_nucleic_acid_sequence(
                    entry.sequence, mol_type.value
                )
            
            if warnings:
                for warning in warnings:
                    print(f"  Warning: Sequence {i+1} ({entry.identifier}): {warning}")
            
            if not is_valid:
                all_valid = False
                print(f"  Error: Sequence {i+1} ({entry.identifier}) is invalid", file=sys.stderr)
        
        if not all_valid:
            print("Error: Some sequences are invalid. Use --no-validate to skip validation.", file=sys.stderr)
            sys.exit(1)
    elif use_mixed_mode and not args.no_validate and args.verbose:
        print("Validation: Skipped (mixed mode validates per-sequence during type detection)")
    
    # Generate name if not provided
    if args.name is None:
        name = generate_name(fasta_entries, args.fasta)
        if args.verbose:
            print(f"Generated target name: {name}")
    else:
        name = args.name
    
    # Create AlphaFold 3 input
    if args.verbose:
        print(f"Creating AlphaFold 3 input with {len(seeds)} seed(s)...")
    
    if use_mixed_mode:
        # Mixed mode: auto-detect type for each sequence
        af3_input = create_af3_input_mixed(
            fasta_entries=fasta_entries,
            name=name,
            seeds=seeds,
            verbose=args.verbose,
            version=args.version,
        )
    else:
        # Single mode: use specified or detected type for all sequences
        af3_input = create_af3_input(
            fasta_entries=fasta_entries,
            name=name,
            mol_type=mol_type,
            seeds=seeds,
            version=args.version,
        )
    
    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(af3_input, f, indent=2)
    
    print(f"Successfully created AlphaFold 3 input: {output_path}")
    print(f"  Target name: {name}")
    print(f"  Sequences: {len(fasta_entries)}")
    print(f"  Seeds: {seeds}")
    print(f"  Total residues: {sum(len(e.sequence) for e in fasta_entries)}")
    
    if use_mixed_mode:
        # Count types in output
        type_counts = {}
        for seq_entry in af3_input['sequences']:
            for mol_type in seq_entry.keys():
                type_counts[mol_type] = type_counts.get(mol_type, 0) + 1
        
        print(f"  Molecule types (auto-detected per sequence):")
        for mol_type, count in sorted(type_counts.items()):
            print(f"    {mol_type}: {count} chain(s)")
    else:
        print(f"  Molecule type: {mol_type.value}")
    
    if len(fasta_entries) == 1:
        print(f"  Chain: A")
    else:
        chains = assign_chains(len(fasta_entries))
        print(f"  Chains: {', '.join(chains)}")


if __name__ == '__main__':
    main()

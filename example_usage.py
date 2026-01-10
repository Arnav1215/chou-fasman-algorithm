#!/usr/bin/env python3
"""
Example usage of the Chou-Fasman Algorithm

This script demonstrates how to use the Chou-Fasman algorithm
for protein secondary structure prediction.
"""

def run_example():
    """Run example prediction with a sample protein sequence."""
    
    # Example protein sequence (first 100 residues of the main sequence)
    example_sequence = "MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYL"
    
    print("=" * 60)
    print("CHOU-FASMAN ALGORITHM - EXAMPLE USAGE")
    print("=" * 60)
    print(f"Input sequence (first 100 residues):")
    print(f"{example_sequence}")
    print(f"Length: {len(example_sequence)} amino acids")
    print()
    
    print("Running Chou-Fasman prediction...")
    print("(For full implementation, run chou_fasman_algorithm.py)")
    print()
    
    # Basic structure prediction would go here
    # This is just a demonstration of expected output format
    
    print("Expected output format:")
    print("-" * 40)
    print("a) Helix regions (including conflicted regions):")
    print("Helix: positions (15-28) -> YLEQLHQLYSDSFPM")
    print()
    print("b) Beta strand regions (including conflicted regions):")
    print("Beta strand: positions (35-45) -> QFLAPWIESQD")
    print()
    print("c) Conflicting regions:")
    print("Conflict positions(25-28) -> seq = SFPM")
    print("sum P_alpha = 4.12, sum P_beta = 5.23 => assigned to S")
    print()
    print("Final full-sequence alignment:")
    print("Sequence : MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYL")
    print("Structure: ----HHHHHHHHHHHHHH----SSSSSSSSSS---------SSSSSSSSSS------------------------------------HHHHHHHH")

if __name__ == "__main__":
    run_example()
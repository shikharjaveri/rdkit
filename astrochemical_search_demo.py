#!/usr/bin/env python
"""
Demonstration of using the modified RDKit for substructure searching
in astrochemical molecules and generating chemical fingerprints.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFingerprintGenerator
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit.Chem.Draw import rdMolDraw2D

def simulate_astro_sanitize(mol):
    """
    Python implementation that simulates our C++ sanitizeAstroMol function
    Uses only the operations that are safe for astrochemical molecules
    """
    # First, make sure we're working with a writable molecule
    if not isinstance(mol, Chem.RWMol):
        mol = Chem.RWMol(mol)
    
    # Update property cache without checking valence - this is key for unusual molecules
    mol.UpdatePropertyCache(strict=False)
    
    # Symmetrize rings
    try:
        Chem.GetSymmSSSR(mol)
    except Exception as e:
        print(f"  Warning: GetSymmSSSR failed - {str(e)}")
    
    # Set aromaticity
    try:
        Chem.SetAromaticity(mol)
    except Exception as e:
        print(f"  Warning: SetAromaticity failed - {str(e)}")
    
    # Set conjugation
    try:
        Chem.SetConjugation(mol)
    except Exception as e:
        print(f"  Warning: SetConjugation failed - {str(e)}")
    
    # Set hybridization
    try:
        Chem.SetHybridization(mol)
    except Exception as e:
        print(f"  Warning: SetHybridization failed - {str(e)}")
    
    # Update property cache again
    mol.UpdatePropertyCache(strict=False)
    
    return mol

def create_astro_mol(smiles):
    """
    Create a molecule using our astrochemical sanitization approach
    """
    if hasattr(Chem, 'MolFromSmilesAstro'):
        return Chem.MolFromSmilesAstro(smiles)
    else:
        # Simulate the function if not available
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        return simulate_astro_sanitize(mol)

def highlight_substructure(mol, smarts_pattern, filename):
    """
    Generate an image of a molecule with highlighted substructure matches
    """
    patt = Chem.MolFromSmarts(smarts_pattern)
    if patt is None:
        print(f"Error: Invalid SMARTS pattern: {smarts_pattern}")
        return
    
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        print(f"No matches found for pattern: {smarts_pattern}")
        return
    
    # Prepare coordinates for visualization
    AllChem.Compute2DCoords(mol)
    
    # Create the drawing
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 300)
    drawer.drawOptions().addAtomIndices = True
    
    # Prepare highlighting
    atoms_to_highlight = set()
    bonds_to_highlight = set()
    colors = {}
    
    for match in matches:
        # Add atoms to highlight
        for atom_idx in match:
            atoms_to_highlight.add(atom_idx)
            colors[atom_idx] = (1, 0, 0)  # Red color for matching atoms
        
        # Add bonds between matched atoms to highlight
        for i in range(len(match)-1):
            for j in range(i+1, len(match)):
                atom1, atom2 = match[i], match[j]
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond is not None:
                    bond_idx = bond.GetIdx()
                    bonds_to_highlight.add(bond_idx)
    
    # Draw with highlighting  
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, 
                                       highlightAtoms=list(atoms_to_highlight),
                                       highlightBonds=list(bonds_to_highlight),
                                       highlightAtomColors=colors)
    drawer.FinishDrawing()
    
    # Save the image
    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())
    
    print(f"Saved highlighted structure to {filename}")
    print(f"Found {len(matches)} matches for pattern: {smarts_pattern}")
    
    return matches

def define_astrochemical_molecules():
    """
    Define a set of astrochemical molecules to work with
    """
    molecules = [
        # Simple molecules
        ("Polyyne anion", "C#CC#CC#CC#CC#[C-]"),
        ("Cyanomethylidyne radical", "CC#N[C]"),
        ("Carbon monoxide", "[C]=O"),
        ("Carbon dioxide", "[C](=O)=O"),
        
        # Astrochemical ions
        ("HCO+", "[CH+]=O"),
        ("N2H+", "[NH+]#N"),
        
        # PAHs and derivatives
        ("Naphthalene", "c1ccc2ccccc2c1"),
        ("Radical naphthalene", "c1ccc2c(ccc[c]2)c1"),
        
        # Larger structures
        ("Pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
        ("Polyynic chain", "C#CC#CC#CC#CC#CC#CC#CC#C")
    ]
    return molecules

def define_smarts_patterns():
    """
    Define SMARTS patterns for substructure searching in astrochemical molecules
    """
    patterns = [
        # Radical patterns
        ("Radical carbon", "[C]"),
        ("Radical with triple bond", "[C]#C"),
        
        # Carbon chain patterns
        ("C≡C-C≡C motif", "C#CC#C"),
        ("Long polyyne segment", "C#CC#CC#C"),
        
        # Ionic patterns
        ("Positively charged nitrogen", "[N+]"),
        ("Positively charged carbon", "[C+]"),
        
        # Hypervalent patterns
        ("Hypervalent carbon", "[C](=*)=*"),
        ("Terminal triple bond", "*#[CH]"),
        
        # PAH patterns
        ("Fused benzene rings", "c1ccc2ccccc2c1")
    ]
    return patterns

def generate_fingerprints(mols, names):
    """
    Generate and compare fingerprints for a set of molecules
    """
    # Create fingerprint generator
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    
    fingerprints = []
    for mol in mols:
        if mol is not None:
            fp = morgan_gen.GetFingerprint(mol)
            fingerprints.append(fp)
        else:
            fingerprints.append(None)
    
    # Calculate Tanimoto similarities
    n = len(fingerprints)
    similarity_matrix = np.zeros((n, n))
    
    for i in range(n):
        if fingerprints[i] is None:
            continue
        for j in range(n):
            if fingerprints[j] is None:
                similarity_matrix[i, j] = 0
            else:
                similarity_matrix[i, j] = rdFingerprintGenerator.GetTanimotoSimilarity(fingerprints[i], fingerprints[j])
    
    # Create and display similarity matrix as a pandas DataFrame
    df = pd.DataFrame(similarity_matrix, index=names, columns=names)
    print("\nTanimoto Similarity Matrix:")
    print(df)
    
    # Save the matrix to a CSV file
    df.to_csv("tanimoto_similarity.csv")
    print("Saved similarity matrix to tanimoto_similarity.csv")
    
    return fingerprints, similarity_matrix

def plot_similarity_heatmap(similarity_matrix, names, filename="similarity_heatmap.png"):
    """
    Generate a heatmap of the similarity matrix
    """
    plt.figure(figsize=(10, 8))
    plt.imshow(similarity_matrix, cmap="YlGnBu", interpolation='nearest')
    plt.colorbar(label='Tanimoto Similarity')
    
    # Add labels
    plt.xticks(range(len(names)), names, rotation=45, ha='right')
    plt.yticks(range(len(names)), names)
    
    # Add similarity values in the cells
    for i in range(len(names)):
        for j in range(len(names)):
            text = plt.text(j, i, f"{similarity_matrix[i, j]:.2f}",
                          ha="center", va="center", color="black" if similarity_matrix[i, j] < 0.7 else "white")
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Saved similarity heatmap to {filename}")

def main():
    """
    Main function to demonstrate searching for substructures in astrochemical molecules
    """
    print("="*80)
    print("RDKIT ASTROCHEMICAL SUBSTRUCTURE SEARCH DEMONSTRATION")
    print("="*80)
    
    # Create molecules
    molecules = define_astrochemical_molecules()
    mols = []
    mol_names = []
    
    for name, smiles in molecules:
        print(f"\nProcessing: {name} ({smiles})")
        mol = create_astro_mol(smiles)
        
        if mol is not None:
            print(f"  ✓ Created molecule successfully")
            print(f"    SMILES: {Chem.MolToSmiles(mol)}")
            
            # Generate 2D coordinates for the molecule
            AllChem.Compute2DCoords(mol)
            
            # Save an image of the molecule
            img_file = f"{name.lower().replace(' ', '_')}.png"
            img = Draw.MolToImage(mol, size=(300, 200))
            img.save(img_file)
            print(f"  ✓ Saved image to {img_file}")
            
            mols.append(mol)
            mol_names.append(name)
        else:
            print(f"  ✗ Failed to create molecule")
            mols.append(None)
            mol_names.append(name)
    
    # Perform substructure searches
    smarts_patterns = define_smarts_patterns()
    
    print("\n" + "="*60)
    print("SUBSTRUCTURE SEARCHING")
    print("="*60)
    
    # Select a subset of molecules to demonstrate substructure searching
    demo_mols = [
        ("Polyyne anion", mols[0]),
        ("Radical naphthalene", mols[7]),
        ("Carbon dioxide", mols[3])
    ]
    
    for pattern_name, smarts in smarts_patterns:
        print(f"\nSearching for pattern: {pattern_name} ({smarts})")
        
        for mol_name, mol in demo_mols:
            if mol is not None:
                print(f"\n  Testing {mol_name}:")
                
                # Highlight and save image with matches
                img_file = f"{mol_name.lower().replace(' ', '_')}_{pattern_name.lower().replace(' ', '_')}.png"
                matches = highlight_substructure(mol, smarts, img_file)
                
                if matches:
                    print(f"    Found {len(matches)} matches")
                    for idx, match in enumerate(matches):
                        print(f"    Match {idx+1}: atoms {match}")
                else:
                    print(f"    No matches found")
    
    # Generate fingerprints and calculate similarities
    print("\n" + "="*60)
    print("FINGERPRINT GENERATION AND SIMILARITY COMPARISON")
    print("="*60)
    
    # Filter out None molecules
    valid_mols = []
    valid_names = []
    for mol, name in zip(mols, mol_names):
        if mol is not None:
            valid_mols.append(mol)
            valid_names.append(name)
    
    fingerprints, similarity_matrix = generate_fingerprints(valid_mols, valid_names)
    
    # Plot similarity heatmap
    plot_similarity_heatmap(similarity_matrix, valid_names)
    
    print("\n" + "="*80)
    print("Demonstration completed successfully!")
    print("="*80)

if __name__ == "__main__":
    main() 
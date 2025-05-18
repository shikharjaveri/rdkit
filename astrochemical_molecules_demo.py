#!/usr/bin/env python
"""
Demonstration of the new RDKit astrochemical molecule sanitization features.
This script shows how to use the custom sanitization workflow on various
challenging molecules that would normally fail with standard RDKit sanitization.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def create_mol_grid(mols, legends, filename='molecule_grid.png', mol_per_row=3):
    """Create a grid of molecule images with legends"""
    if len(mols) == 0:
        return
    
    n_mols = len(mols)
    n_rows = (n_mols + mol_per_row - 1) // mol_per_row
    
    # Set up figure
    fig = plt.figure(figsize=(4*mol_per_row, 3.5*n_rows))
    gs = GridSpec(n_rows, mol_per_row)
    
    for i, (mol, legend) in enumerate(zip(mols, legends)):
        row = i // mol_per_row
        col = i % mol_per_row
        
        img = Draw.MolToImage(mol, size=(300, 200))
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(img)
        ax.set_title(legend, fontsize=12)
        ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=100, bbox_inches='tight')
    print(f"Molecule grid saved as {filename}")

def load_molecules():
    """Define a set of astrochemical molecules that are challenging for standard RDKit"""
    molecules = [
        # Astrochemical molecules with unusual valence
        ("Polyyne anion", "C#CC#CC#CC#CC#[C-]"),
        ("Cyanide radical", "[C]#N"),
        ("Carbon monoxide with unusual valence", "[C]=O"),
        ("Cyanomethylidyne radical", "CC#N[C]"),
        
        # Ions common in astrochemistry
        ("HCO+ (formyl cation)", "[CH+]=O"),
        ("N2H+ (diazenylium)", "[NH+]#N"),
        ("H3+ (trihydrogen cation)", "[H+]([H])[H]"),
        
        # Other unusual structures
        ("Carbon dioxide (unusual representation)", "[C](=O)(=O)"),
        ("Silicon dioxide", "O=[Si]=O"),
        ("Hypervalent sulfur species", "O=S(=O)(=O)[O-]"),
        
        # PAHs and derivatives
        ("Benzyne (strained triple bond)", "C1=CC=CC=C#1"),
        ("Radical PAH fragment", "c1cc[c]cc1")
    ]
    return molecules

def process_with_standard_rdkit(smiles, name):
    """Try to process a molecule with standard RDKit sanitization"""
    print(f"\nStandard RDKit approach for {name} ({smiles}):")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            print("  ✓ Molecule created successfully with standard sanitization")
            print(f"    SMILES: {Chem.MolToSmiles(mol)}")
            return mol
        else:
            print("  ✗ Failed to create molecule with standard sanitization")
            return None
    except Exception as e:
        print(f"  ✗ Error during standard sanitization: {e}")
        return None

def process_with_workaround(smiles, name):
    """Process a molecule with the workaround approach (no sanitization + property cache update)"""
    print(f"\nWorkaround approach for {name} ({smiles}):")
    try:
        # Create the molecule without sanitization
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            print("  ✗ Failed to create initial molecule structure")
            return None
        
        # Apply the workaround approach
        print("  → Updating property cache with strict=False")
        mol.UpdatePropertyCache(strict=False)
        
        print("  → Applying safe sanitization operations")
        Chem.GetSymmSSSR(mol)
        Chem.SetAromaticity(mol)
        Chem.SetConjugation(mol)
        Chem.SetHybridization(mol)
        mol.UpdatePropertyCache(strict=False)
        
        print("  ✓ Molecule processed successfully with workaround")
        print(f"    SMILES: {Chem.MolToSmiles(mol)}")
        return mol
    except Exception as e:
        print(f"  ✗ Error during workaround approach: {e}")
        return None

def process_with_astrochemical_function(smiles, name):
    """Process a molecule with our new astrochemical sanitization function"""
    print(f"\nNew astrochemical approach for {name} ({smiles}):")
    try:
        # Use our new function if available, otherwise simulate it
        if hasattr(Chem, 'MolFromSmilesAstro'):
            print("  → Using MolFromSmilesAstro")
            mol = Chem.MolFromSmilesAstro(smiles)
        else:
            print("  → Simulating sanitizeAstroMol function")
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            mol.UpdatePropertyCache(strict=False)
            
            # Apply only safe sanitization operations
            Chem.GetSymmSSSR(mol)
            Chem.SetAromaticity(mol)
            Chem.SetConjugation(mol)
            Chem.SetHybridization(mol)
            mol.UpdatePropertyCache(strict=False)
        
        if mol is None:
            print("  ✗ Failed to create molecule")
            return None
        
        print("  ✓ Molecule processed successfully with astrochemical sanitization")
        print(f"    SMILES: {Chem.MolToSmiles(mol)}")
        return mol
    except Exception as e:
        print(f"  ✗ Error during astrochemical sanitization: {e}")
        return None

def generate_3d_structure(mol, name):
    """Try to generate a 3D structure for the molecule"""
    if mol is None:
        return None
    
    try:
        # Make a copy to avoid modifying the original
        mol_3d = Chem.Mol(mol)
        AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol_3d)
        print(f"  ✓ Generated 3D structure for {name}")
        return mol_3d
    except Exception as e:
        print(f"  ✗ Failed to generate 3D structure for {name}: {e}")
        return None

def calculate_properties(mol, name):
    """Calculate and display properties for the molecule"""
    if mol is None:
        return
    
    print(f"\nProperties for {name}:")
    try:
        # Basic properties
        print(f"  Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
        print(f"  Atoms: {mol.GetNumAtoms()}")
        print(f"  Bonds: {mol.GetNumBonds()}")
        
        # Calculate formal charge
        formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        print(f"  Formal charge: {formal_charge}")
        
        # Calculate radical electrons
        radical_e = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        print(f"  Radical electrons: {radical_e}")
        
        # Atom-specific properties
        print("\n  Atom properties:")
        for atom in mol.GetAtoms():
            print(f"    Atom {atom.GetIdx()} ({atom.GetSymbol()}):")
            print(f"      Hybridization: {atom.GetHybridization()}")
            print(f"      Formal charge: {atom.GetFormalCharge()}")
            print(f"      Explicit valence: {atom.GetExplicitValence()}")
            print(f"      Implicit valence: {atom.GetImplicitValence()}")
    except Exception as e:
        print(f"  ✗ Error calculating properties: {e}")

def main():
    """Main function to demonstrate using RDKit with astrochemical molecules"""
    print("="*80)
    print("RDKIT ASTROCHEMICAL MOLECULES DEMONSTRATION")
    print("="*80)
    
    molecules = load_molecules()
    
    # Lists to store successfully processed molecules
    standard_mols = []
    astro_mols = []
    standard_legends = []
    astro_legends = []
    
    # Process each molecule
    for name, smiles in molecules:
        print("\n" + "="*60)
        print(f"PROCESSING: {name} ({smiles})")
        print("="*60)
        
        # Try standard RDKit approach
        std_mol = process_with_standard_rdkit(smiles, name)
        if std_mol:
            std_mol_2d = std_mol
            AllChem.Compute2DCoords(std_mol_2d)
            standard_mols.append(std_mol_2d)
            standard_legends.append(f"{name} (standard)")
        
        # Try astrochemical approach
        astro_mol = process_with_astrochemical_function(smiles, name)
        if astro_mol:
            astro_mol_2d = astro_mol
            AllChem.Compute2DCoords(astro_mol_2d)
            astro_mols.append(astro_mol_2d)
            astro_legends.append(f"{name} (astro)")
            
            # Calculate properties for the astrochemical molecule
            calculate_properties(astro_mol, name)
            
            # Try to generate 3D structure
            mol_3d = generate_3d_structure(astro_mol, name)
            if mol_3d:
                # Save as MOL file if successful
                mol_file = f"{name.lower().replace(' ', '_')}_3d.mol"
                Chem.MolToMolFile(mol_3d, mol_file)
                print(f"  ✓ Saved 3D structure to {mol_file}")
    
    # Generate images with all the molecules
    if standard_mols:
        create_mol_grid(standard_mols, standard_legends, "standard_molecules.png")
    if astro_mols:
        create_mol_grid(astro_mols, astro_legends, "astrochemical_molecules.png")
    
    print("\n" + "="*80)
    print(f"Successfully processed {len(astro_mols)} molecules with astrochemical approach")
    print(f"Only {len(standard_mols)} molecules succeeded with standard RDKit sanitization")
    print("="*80)

if __name__ == "__main__":
    main() 
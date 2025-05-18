#!/usr/bin/env python
"""
A simplified script to show properties of astrochemical molecules after custom sanitization
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

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

def print_property(name, value, indent=0):
    """Helper function to nicely print property name and value"""
    spaces = ' ' * indent
    print(f"{spaces}{name}: {value}")

def get_molecule_properties(mol):
    """Print the main properties of a molecule"""
    if mol is None:
        print("Invalid molecule - cannot calculate properties")
        return
    
    print("\n" + "="*70)
    print("PROPERTIES ACCESSIBLE AFTER ASTROCHEMICAL SANITIZATION")
    print("="*70)
    
    # BASIC PROPERTIES
    print("\n1. BASIC MOLECULAR PROPERTIES")
    print_property("SMILES", Chem.MolToSmiles(mol), 2)
    print_property("Formula", Chem.rdMolDescriptors.CalcMolFormula(mol), 2)
    print_property("Number of atoms", mol.GetNumAtoms(), 2)
    print_property("Number of heavy atoms", Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol), 2)
    print_property("Number of bonds", mol.GetNumBonds(), 2)
    
    # Calculate formal charge as sum of all atom formal charges
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    print_property("Formal charge", formal_charge, 2)
    
    # Calculate radical electrons as sum from all atoms
    rad_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    print_property("Number of radical electrons", rad_electrons, 2)
    
    # ATOMIC PROPERTIES
    print("\n2. ATOMIC PROPERTIES")
    for atom in mol.GetAtoms():
        print(f"  Atom {atom.GetIdx()} ({atom.GetSymbol()}):")
        print_property("Atomic Number", atom.GetAtomicNum(), 4)
        print_property("Formal Charge", atom.GetFormalCharge(), 4)
        print_property("Hybridization", str(atom.GetHybridization()), 4)
        print_property("Aromatic", atom.GetIsAromatic(), 4)
        print_property("In Ring", atom.IsInRing(), 4)
        print_property("Explicit Valence", atom.GetExplicitValence(), 4)
        print_property("Implicit Valence", atom.GetImplicitValence(), 4)
        print_property("Total Valence", atom.GetTotalValence(), 4)
        print_property("Num Explicit Hs", atom.GetNumExplicitHs(), 4)
        print_property("Num Implicit Hs", atom.GetNumImplicitHs(), 4)
        print_property("Num Radical Electrons", atom.GetNumRadicalElectrons(), 4)
    
    # BOND PROPERTIES
    print("\n3. BOND PROPERTIES")
    for bond in mol.GetBonds():
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        print(f"  Bond {bond.GetIdx()} ({begin_atom.GetSymbol()}-{end_atom.GetSymbol()}):")
        print_property("Bond Type", str(bond.GetBondType()), 4)
        print_property("Is Conjugated", bond.GetIsConjugated(), 4)
        print_property("Is Aromatic", bond.GetIsAromatic(), 4)
        print_property("Is In Ring", bond.IsInRing(), 4)
    
    print("\n" + "="*70)
    
    # Try to generate and save molecular structure
    try:
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(400, 200))
        img.save('molecule_structure.png')
        print("Saved molecular structure as 'molecule_structure.png'")
    except Exception as e:
        print(f"Could not generate 2D structure: {str(e)}")

def main():
    smiles = "C#CC#CC#CC#CC#[C-]"  # Polyyne anion
    print(f"\nTesting molecule: {smiles}")
    
    # Create molecule without sanitization
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    
    if mol:
        # Apply our custom sanitization
        mol = simulate_astro_sanitize(mol)
        
        # Print properties
        get_molecule_properties(mol)
    else:
        print(f"Failed to create molecule for {smiles}")

if __name__ == "__main__":
    main() 
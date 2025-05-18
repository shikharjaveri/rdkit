#!/usr/bin/env python
"""
Show properties of an astrochemical molecule using our custom sanitization approach
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

def print_property(name, value):
    """Helper function to nicely print property name and value"""
    print(f"  {name}: {value}")

def show_molecule_properties(mol):
    """Print the main properties of a molecule"""
    if mol is None:
        print("Invalid molecule - cannot calculate properties")
        return
    
    print("\nBASIC PROPERTIES:")
    print_property("SMILES", Chem.MolToSmiles(mol))
    print_property("Formula", Chem.rdMolDescriptors.CalcMolFormula(mol))
    print_property("Number of atoms", mol.GetNumAtoms())
    print_property("Number of heavy atoms", Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol))
    print_property("Number of bonds", mol.GetNumBonds())
    
    # Calculate formal charge as sum of all atom formal charges
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    print_property("Formal charge", formal_charge)
    
    # Calculate radical electrons as sum from all atoms
    rad_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    print_property("Number of radical electrons", rad_electrons)
    
    print("\nATOMIC PROPERTIES:")
    for atom in mol.GetAtoms():
        print(f"  Atom {atom.GetIdx()} ({atom.GetSymbol()}):")
        print(f"    Atomic Number: {atom.GetAtomicNum()}")
        print(f"    Formal Charge: {atom.GetFormalCharge()}")
        print(f"    Hybridization: {atom.GetHybridization()}")
        print(f"    Explicit Valence: {atom.GetExplicitValence()}")
        print(f"    Implicit Valence: {atom.GetImplicitValence()}")
        print(f"    Total Valence: {atom.GetTotalValence()}")
        print(f"    Num Explicit Hs: {atom.GetNumExplicitHs()}")
        print(f"    Num Implicit Hs: {atom.GetNumImplicitHs()}")
        print(f"    Num Radical Electrons: {atom.GetNumRadicalElectrons()}")
    
    print("\nBOND PROPERTIES:")
    for bond in mol.GetBonds():
        begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        print(f"  Bond {bond.GetIdx()} ({begin_atom.GetSymbol()}-{end_atom.GetSymbol()}):")
        print(f"    Bond Type: {bond.GetBondType()}")
        print(f"    Is Conjugated: {bond.GetIsConjugated()}")
        print(f"    Is Aromatic: {bond.GetIsAromatic()}")
        print(f"    Is In Ring: {bond.IsInRing()}")

def main():
    # The polyyne anion
    smiles = "C#CC#CC#CC#CC#[C-]"
    print(f"TESTING MOLECULE: {smiles}")
    
    # Create molecule without sanitization
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    
    if mol:
        # Apply our custom sanitization
        mol = simulate_astro_sanitize(mol)
        
        # Print properties
        show_molecule_properties(mol)
        
        # Try to generate 2D coordinates for visualization
        try:
            AllChem.Compute2DCoords(mol)
            img = Draw.MolToImage(mol, size=(400, 200))
            img.save('polyyne_structure.png')
            print("\nSaved molecular structure as 'polyyne_structure.png'")
        except Exception as e:
            print(f"\nCould not generate 2D structure: {str(e)}")
    else:
        print(f"Failed to create molecule for {smiles}")

if __name__ == "__main__":
    main() 
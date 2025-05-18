#!/usr/bin/env python
"""
Test script for astrochemical molecule sanitization
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import sys

def test_astro_mol_sanitization():
    """Test that astrochemical molecules can be sanitized"""
    # Examples of astrochemical molecules with unusual valence states
    # that would normally fail sanitization
    smiles_list = [
        'CC#N[C]',          # Cyanomethylidyne radical
        'C#C[C]',           # Propargyl radical
        'N#N#N',            # Azide with unusual bonding
        '[C](=O)(=O)',      # Carbon dioxide with unusual representation
        'O=[Si]=O',         # Silicon dioxide
        'CC(=O)[O]',        # Acetate radical
    ]
    
    success_count = 0
    
    for idx, smiles in enumerate(smiles_list):
        print(f"\nTesting molecule {idx+1}: {smiles}")
        
        # First show that standard sanitization fails
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            Chem.SanitizeMol(mol)
            print("  Standard sanitization: SUCCEEDED (unexpected!)")
        except:
            print("  Standard sanitization: FAILED (expected)")
        
        # Now try with our custom sanitization
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            Chem.SanitizeAstroMol(mol)
            success_count += 1
            print("  Astro sanitization: SUCCEEDED")
            
            # Test that we can now use other RDKit functionality
            print("  Molecule has properties:")
            print(f"    - Num atoms: {mol.GetNumAtoms()}")
            print(f"    - Num bonds: {mol.GetNumBonds()}")
            print(f"    - Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
            
            # Try to generate 3D coordinates
            try:
                status = AllChem.EmbedMolecule(mol)
                if status == 0:
                    print("  3D coordinate generation: SUCCEEDED")
                else:
                    print("  3D coordinate generation: FAILED")
            except Exception as e:
                print(f"  3D coordinate generation: ERROR - {str(e)}")
            
        except Exception as e:
            print(f"  Astro sanitization: FAILED - {str(e)}")
    
    print(f"\nSummary: {success_count}/{len(smiles_list)} molecules passed astro sanitization")
    
    return success_count == len(smiles_list)

if __name__ == "__main__":
    success = test_astro_mol_sanitization()
    sys.exit(0 if success else 1) 
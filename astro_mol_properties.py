#!/usr/bin/env python
"""
Demonstrates all properties that can be obtained from an astrochemical molecule
after applying our custom sanitization approach.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, Lipinski, rdMolDescriptors, QED
import numpy as np
import sys

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

def get_all_properties(mol):
    """
    Calculate and print all available properties for an astrochemical molecule
    that has been sanitized using our custom approach
    """
    if mol is None:
        print("Invalid molecule - cannot calculate properties")
        return
    
    print("\n" + "="*70)
    print("PROPERTIES OF ASTROCHEMICAL MOLECULE AFTER CUSTOM SANITIZATION")
    print("="*70)
    
    # BASIC PROPERTIES
    print("\n1. BASIC MOLECULAR PROPERTIES")
    print_property("SMILES", Chem.MolToSmiles(mol), 2)
    print_property("Formula", Chem.rdMolDescriptors.CalcMolFormula(mol), 2)
    print_property("Number of atoms", mol.GetNumAtoms(), 2)
    print_property("Number of heavy atoms", Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol), 2)
    print_property("Number of bonds", mol.GetNumBonds(), 2)
    print_property("Number of rotatable bonds", Chem.rdMolDescriptors.CalcNumRotatableBonds(mol), 2)
    
    # Calculate formal charge as sum of all atom formal charges
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    print_property("Formal charge", formal_charge, 2)
    
    try:
        num_radical = Chem.rdMolDescriptors.CalcNumRadicalElectrons(mol)
        print_property("Number of radical electrons", num_radical, 2)
    except:
        rad_electrons = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        print_property("Number of radical electrons", rad_electrons, 2)
        
    print_property("Is chiral", mol.HasProp("_ChiralityPossible"), 2)
    
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
        print_property("Chiral Tag", str(atom.GetChiralTag()), 4)
        print()
    
    # BOND PROPERTIES
    print("\n3. BOND PROPERTIES")
    for bond in mol.GetBonds():
        print(f"  Bond {bond.GetIdx()} ({mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()}-{mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()}):")
        print_property("Bond Type", str(bond.GetBondType()), 4)
        print_property("Is Conjugated", bond.GetIsConjugated(), 4)
        print_property("Is Aromatic", bond.GetIsAromatic(), 4)
        print_property("Is In Ring", bond.IsInRing(), 4)
        print_property("Stereo", str(bond.GetStereo()), 4)
        print()
    
    # PHYSICAL PROPERTIES
    print("\n4. PHYSICAL AND CHEMICAL PROPERTIES")
    try:
        print_property("Molecular Weight", Descriptors.MolWt(mol), 2)
    except Exception as e:
        print_property("Molecular Weight", f"Error: {str(e)}", 2)
    
    try:
        print_property("Exact Mass", Descriptors.ExactMolWt(mol), 2)
    except Exception as e:
        print_property("Exact Mass", f"Error: {str(e)}", 2)
    
    try:
        print_property("Heavy Atom Molecular Weight", Descriptors.HeavyAtomMolWt(mol), 2)
    except Exception as e:
        print_property("Heavy Atom Molecular Weight", f"Error: {str(e)}", 2)
    
    try:
        print_property("LogP", Descriptors.MolLogP(mol), 2)
    except Exception as e:
        print_property("LogP", f"Error: {str(e)}", 2)
    
    try:
        print_property("TPSA", Descriptors.TPSA(mol), 2)
    except Exception as e:
        print_property("TPSA", f"Error: {str(e)}", 2)
    
    try:
        print_property("LabuteASA", Descriptors.LabuteASA(mol), 2)
    except Exception as e:
        print_property("LabuteASA", f"Error: {str(e)}", 2)
    
    try:
        print_property("Num Hydrogen Bond Donors", Descriptors.NumHDonors(mol), 2)
    except Exception as e:
        print_property("Num Hydrogen Bond Donors", f"Error: {str(e)}", 2)
    
    try:
        print_property("Num Hydrogen Bond Acceptors", Descriptors.NumHAcceptors(mol), 2)
    except Exception as e:
        print_property("Num Hydrogen Bond Acceptors", f"Error: {str(e)}", 2)
    
    # TOPOLOGICAL PROPERTIES
    print("\n5. TOPOLOGICAL PROPERTIES")
    try:
        print_property("Number of Rings", Chem.rdMolDescriptors.CalcNumRings(mol), 2)
    except Exception as e:
        print_property("Number of Rings", f"Error: {str(e)}", 2)
    
    try:
        print_property("Number of Aromatic Rings", Chem.rdMolDescriptors.CalcNumAromaticRings(mol), 2)
    except Exception as e:
        print_property("Number of Aromatic Rings", f"Error: {str(e)}", 2)
    
    try:
        print_property("Number of Aliphatic Rings", Chem.rdMolDescriptors.CalcNumAliphaticRings(mol), 2)
    except Exception as e:
        print_property("Number of Aliphatic Rings", f"Error: {str(e)}", 2)
    
    try:
        print_property("Number of Saturated Rings", Chem.rdMolDescriptors.CalcNumSaturatedRings(mol), 2)
    except Exception as e:
        print_property("Number of Saturated Rings", f"Error: {str(e)}", 2)
    
    try:
        print_property("Number of Heterocycles", Chem.rdMolDescriptors.CalcNumHeterocycles(mol), 2)
    except Exception as e:
        print_property("Number of Heterocycles", f"Error: {str(e)}", 2)
    
    try:
        print_property("Chi0v", Descriptors.Chi0v(mol), 2)
    except Exception as e:
        print_property("Chi0v", f"Error: {str(e)}", 2)
    
    try:
        print_property("Chi1v", Descriptors.Chi1v(mol), 2)
    except Exception as e:
        print_property("Chi1v", f"Error: {str(e)}", 2)
    
    try:
        print_property("Kappa1", Descriptors.Kappa1(mol), 2)
    except Exception as e:
        print_property("Kappa1", f"Error: {str(e)}", 2)
    
    try:
        print_property("Kappa2", Descriptors.Kappa2(mol), 2)
    except Exception as e:
        print_property("Kappa2", f"Error: {str(e)}", 2)
    
    # ELECTRONIC PROPERTIES
    print("\n6. ELECTRONIC PROPERTIES")
    try:
        print_property("Number of Valence Electrons", Chem.rdMolDescriptors.CalcNumValenceElectrons(mol), 2)
    except Exception as e:
        print_property("Number of Valence Electrons", f"Error: {str(e)}", 2)
    
    try:
        print_property("Fraction of sp3 Carbon Atoms", Chem.rdMolDescriptors.CalcFractionCSP3(mol), 2)
    except Exception as e:
        print_property("Fraction of sp3 Carbon Atoms", f"Error: {str(e)}", 2)
    
    try:
        print_property("PEOE_VSA1", Descriptors.PEOE_VSA1(mol), 2)
    except Exception as e:
        print_property("PEOE_VSA1", f"Error: {str(e)}", 2)
    
    try:
        print_property("SMR_VSA1", Descriptors.SMR_VSA1(mol), 2)
    except Exception as e:
        print_property("SMR_VSA1", f"Error: {str(e)}", 2)
    
    try:
        print_property("SlogP_VSA1", Descriptors.SlogP_VSA1(mol), 2)
    except Exception as e:
        print_property("SlogP_VSA1", f"Error: {str(e)}", 2)
    
    # FINGERPRINTS
    print("\n7. FINGERPRINTS (sample bits shown)")
    try:
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
        print_property("Morgan Fingerprint (ECFP4)", str(list(fp.GetOnBits())[:10]) + "...", 2)
    except Exception as e:
        print_property("Morgan Fingerprint (ECFP4)", f"Error: {str(e)}", 2)
    
    try:
        fp = Chem.RDKFingerprint(mol)
        print_property("RDKit Fingerprint", str(list(fp.GetOnBits())[:10]) + "...", 2)
    except Exception as e:
        print_property("RDKit Fingerprint", f"Error: {str(e)}", 2)
    
    try:
        fp = Chem.rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol)
        print_property("Atom Pair Fingerprint", str(list(fp.GetOnBits())[:10]) + "...", 2)
    except Exception as e:
        print_property("Atom Pair Fingerprint", f"Error: {str(e)}", 2)
    
    # DRUG-LIKE PROPERTIES
    print("\n8. DRUG-LIKE PROPERTIES (if applicable to your use case)")
    try:
        print_property("Lipinski Rule of 5 Failures", Lipinski.NumRotatableBonds(mol), 2)
    except Exception as e:
        print_property("Lipinski Rule of 5 Failures", f"Error: {str(e)}", 2)
    
    try:
        print_property("QED (drug-likeness)", QED.default(mol), 2)
    except Exception as e:
        print_property("QED (drug-likeness)", f"Error: {str(e)}", 2)
    
    print("\n" + "="*70)

def main():
    # Example molecules to test
    molecules = [
        ("Polyyne anion", "C#CC#CC#CC#CC#[C-]"),
        ("Cyanomethylidyne radical", "CC#N[C]"),
        ("Carbon dioxide (unusual representation)", "[C](=O)(=O)"),
        ("Silicon dioxide", "O=[Si]=O")
    ]
    
    for name, smiles in molecules:
        print(f"\n\nTESTING MOLECULE: {name} ({smiles})")
        
        # Create molecule without sanitization
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        
        if mol:
            # Apply our custom sanitization
            mol = simulate_astro_sanitize(mol)
            
            # Calculate and print all properties
            get_all_properties(mol)
            
            # Try to generate 2D coordinates for visualization
            try:
                AllChem.Compute2DCoords(mol)
                img = Draw.MolToImage(mol, size=(400, 200), legend=name)
                img_file = name.lower().replace(" ", "_") + ".png"
                img.save(img_file)
                print(f"Saved structure image as '{img_file}'")
            except Exception as e:
                print(f"2D coordinate generation failed: {str(e)}")
        else:
            print(f"Failed to create molecule for {smiles}")

if __name__ == "__main__":
    main() 
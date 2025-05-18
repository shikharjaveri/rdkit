from rdkit import Chem

# Create the polyyne anion molecule without sanitization
smiles = "C#CC#CC#CC#CC#[C-]"
mol = Chem.MolFromSmiles(smiles, sanitize=False)

# Apply custom sanitization
mol.UpdatePropertyCache(strict=False)

# Print basic properties
print(f"SMILES: {Chem.MolToSmiles(mol)}")
print(f"Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"Number of atoms: {mol.GetNumAtoms()}")
print(f"Number of bonds: {mol.GetNumBonds()}")

# Calculate and print formal charge
total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
print(f"Total formal charge: {total_charge}")

# List atom properties
print("\nATOM PROPERTIES:")
for atom in mol.GetAtoms():
    print(f"Atom {atom.GetIdx()} ({atom.GetSymbol()}):")
    print(f"  Formal Charge: {atom.GetFormalCharge()}")
    print(f"  Explicit Valence: {atom.GetExplicitValence()}")
    print(f"  Implicit Valence: {atom.GetImplicitValence()}")
    print(f"  Total Valence: {atom.GetTotalValence()}")

# List bond properties
print("\nBOND PROPERTIES:")
for bond in mol.GetBonds():
    begin_atom = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
    end_atom = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
    print(f"Bond {bond.GetIdx()} ({begin_atom.GetSymbol()}-{end_atom.GetSymbol()}):")
    print(f"  Bond Type: {bond.GetBondType()}") 
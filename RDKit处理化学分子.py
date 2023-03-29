# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 21:19:12 2018

@author: zwn
"""

# conda create -n drug python=3
# conda activate drug
# conda install -c conda-forge rdkit

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


mol = Chem.MolFromSmiles('CNC(=O)C1=CC=CC=C1SC2=CC3=C(C=C2)C(=N[NH]3)\C=C\C4=CC=CC=N4')
fig = Draw.MolToImage(mol, size=(500, 500))
fig.save('./mols.png')

mol = Chem.MolFromSmiles('[nH]1cnc2cncnc21')
Draw.MolToImage(mol, size=(100, 100), kekulize=True,
                wedgeBonds=True, fitImage=False)


mol = Chem.MolFromSmiles('C1=C[N]C=C1')
Draw.MolToImage(mol, size=(100, 100), kekulize=True,
                wedgeBonds=True, fitImage=False)

Draw.ShowMol(mol, size=(300, 300), kekulize=True, 
             wedgeBonds=True, title='RDKit Molecule')

Draw.MolToMPL(mol, size=(100, 100), kekulize=True, 
             wedgeBonds=True)

# Draw.calcAtomGaussians(mol, a=0.03, step=0.02, weights=None)
# fig.axes[0].imshow(z,cmap=cm.gray,
#         interpolation='bilinear',origin='lower',extent=(0,1,0,1))

# fig=Draw.MolToMPL(m)
# contribs=Crippen.rdMolDescriptors._CalcCrippenContribs(m)
# x,y,z=Draw.calcAtomGaussians(m,0.03,step=0.01,weights=logps)
# fig.axes[0].imshow(z,cmap=cm.jet,interpolation='bilinear',
#         origin='lower',extent=(0,1,0,1))


mol_list = Chem.SDMolSupplier('FDA-approved-Drug-kinase-inhibitor75.sdf')
fig = Draw.MolsToGridImage(
    mol_list, 
    molsPerRow=4, 
    subImgSize=(500, 500)
)
fig.save('./mols.png')


# mols = Chem.SDMolSupplier('.sdf/data/drugbank.sdf')
# Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300))

# aspirin
from rdkit import Chem
from rdkit.Chem import Draw
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
Draw.MolToImage(mol, size=(500, 500))

# output mols as sdf format
# f = Chem.SDWriter('./sdf/test_feed_dict.sdf')
# for mol in list(mols):
#     f.write(mol)
# f.close()

















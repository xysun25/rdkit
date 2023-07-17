# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 15:45:00 2021

@author: Administrator
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs
from rdkit.Chem import Draw


try:
    os.mkdir('excel')
except FileExistsError:
    pass


try:
    os.mkdir('figures')
except FileExistsError:
    pass


# 导入数据 smiles 化学分子数据
df = pd.read_excel('./excel/class.xlsx',header=0,sheet_name="drug")


mols = [Chem.MolFromSmiles(i) for i in df["SMILES"]]

mols = [Chem.AddHs(i) for i in mols]

ECFP = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
ECFP_array = np.asarray(ECFP, dtype=float)
ECFP_array.shape



MACCS = [MACCSkeys.GenMACCSKeys(mol) for mol in mols]
MACCS_array = np.asarray(MACCS, dtype=float)
MACCS_array.shape


X = np.concatenate([ECFP_array, MACCS_array], axis=1)

X.shape
# 1024 + 167

#1024 是1024维度的特征向量 # 1x1024的矩阵
#1024 是1024个维度的特征向量
# 167代表的是167个特征(每个特征可能代表的是化学基团是否存在，存在就是1，不存在就是0)

feature = pd.DataFrame(X)
molecule = ['molecule' + str(i).zfill(4) for i in range(1,81)]

name = pd.DataFrame(molecule)

final = pd.concat([name, feature],axis=1)
final.to_excel("./excel/molecule_ECFP.xlsx",index=False)


try:
    os.mkdir('figures')
except FileExistsError:
    pass

img = ECFP_array[0,:].reshape(32,32)
plt.subplots(figsize=(10, 10))
plt.imshow(img)
plt.axis('off')
plt.savefig("./figures/ECFP.jpg", dpi=500)
plt.show()





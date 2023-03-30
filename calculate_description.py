
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

try:
    os.mkdir('model')
except FileExistsError:
    pass

try:
    os.mkdir('figures')
except FileExistsError:
    pass


# 描述符计算模块
# 1.rdkit.Chem.Lipinski模块
# rdkit中提供了许多描述符的计算方法，可用于分子筛选、成药性评估等。
# 以lipinski类药的相关规则为例，
# 可以通过rdkit.Chem.Lipinski模块进行计算。常用的一些性质举例：
# 氢键受体数 NumHAcceptors
# 氢键供体数 NumHDonors
# 可旋转键数 NumRotatableBonds
# 脂肪环数量 NumAliphaticRings
# 芳香环数量 NumAromaticRings
# SP3杂化碳原子比例FractionCSP3
# ……


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Lipinski
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
fig = Draw.MolToImage(mol, size=(500, 500))
fig.save('./mols.png')

Lipinski.NumHAcceptors(mol)
Lipinski.NumHDonors(mol)

# 将描述符名称(key)和值(val)添加到分子属性中：m.SetProp(key, val)
# 获取分子属性：m.GetProp(key)

Ha = Lipinski.NumHAcceptors(mol)
mol.SetProp('Ha', '%s'%Ha)
mol.GetProp('Ha')


# 2.rdkit.Chem.Descriptors 模块
# 大部分的描述符都可以通过rdkit.Chem.Descriptors模块进行计算。
# 该模块也包含了Lipinski的描述符。常用的一些描述符举例：
# 分子量MolWt
# 脂水分配系数MolLogP
# 拓扑极表面积TPSA
# 以计算TPSA为例：Descriptors.TPSA()

from rdkit.Chem import Descriptors

print(Descriptors.TPSA(mol), Descriptors.MolLogP(mol), Descriptors.MolWt(mol))

from rdkit.ML.Descriptors import MoleculeDescriptors
des_list = ['MolWt', 'NumHAcceptors', 'NumHDonors', 'MolLogP', 'NumRotatableBonds']
calculator = MoleculeDescriptors.MolecularDescriptorCalculator(des_list)
calculator.GetDescriptorSummaries()
tuple1 = calculator.CalcDescriptors(mol)


# 计算所有描述符
from rdkit.Chem import Descriptors

var_list = [x[0] for x in Descriptors._descList]

calculator = MoleculeDescriptors.MolecularDescriptorCalculator(var_list)
column_list = calculator.GetDescriptorSummaries()

# 字典
dic = dict(zip(var_list,column_list))
tuple1 = calculator.CalcDescriptors(mol)
value_list = list(tuple1)
len(value_list)
len(column_list)



import pickle
calculator.SaveState('model/descriptor_calculator')

with open('model/descriptor_calculator', 'rb') as f:
    calc = pickle.load(f)

calc.CalcDescriptors(mol)


# 二、原子描述符可视化
# 可以用相似性地图来查看每个原子对描述符的贡献，更多相似性地图的应用可以查看这篇文章的第二部分。
# 1.原子partial charge可视化
# 计算partial charge要复杂一点。计算出的partial charge存储在每个原子的属性中，可以通过GetDoubleProp（浮点数）或GetProp（字符串）来获取。
# partial charge可以表示电子的分布。分子中的化学键是由分布在相连原子周围的电子对组成的。但由于原子的电负性不同，成键电子并不是均匀分布。电负性大的原子吸电子能力强，成键的电子对会更偏向该原子，导致该原子带有负电的性质。相对应的，电负性小的原子则带有正电性质。Partial charge衡量了电子偏向的程度。
# 先计算partial charge：AllChem.ComputeGasteigerCharges()
# 获取某一个原子：GetAtomWithIdx()
# 获取原子的partial charge：GetDoubleProp('_GasteigerCharge')


from rdkit.Chem import AllChem
mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
Draw.MolToImage(mol, size=(500, 500))
AllChem.ComputeGasteigerCharges(mol)
atom = mol.GetAtomWithIdx(0)
atom.GetDoubleProp('_GasteigerCharge')


# 获取每个原子的partial charge，放到contribs中
contribs_list = [round(mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge'), 2) for i in range(mol.GetNumAtoms())]
print(contribs_list)
len(contribs_list)

# 根据给定的权重，生成分子权重图：
# GetSimilarityMapFromWeights(mol, weights, colorMap, contourLines, ...) 
# mol：要绘制的mol对象 weights：权重 colorMap：matplotlib中的色系 contourLines：等高线数量

from rdkit.Chem.Draw import SimilarityMaps
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs_list, contourLines=10)
fig.savefig('figures/descriptor_visulize1.png', dpi=300, bbox_inches='tight')
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs_list, scale=10)
fig.savefig('figures/descriptor_visulize2.png', dpi=300, bbox_inches='tight')

# 2.原子logP可视化
# logP表示脂水分配系数，该值认为与细胞通透性有一定相关性。
# rdkit中提供的Descriptors.MolLogP()方法可以粗略计算logP值，
# 该方法首先做了一个原子分类系统，根据原子及其相连原子的不同而进行分类，
# 再对化学性质相似、logP贡献相似的类别做合并，
# 最终得到了68种精确的原子类别和4种通配类别，并用SMARTS表示。
# 计算时，对一个分子中所有原子进行分类，再乘以每一类的权重并加和，最终得到LogP值。
# 该方法在9920个分子的训练集上的r2为0.918，标准差为0.667。
# 此外摩尔折射率（molar refractivity，MR）也可以通过这种方法计算得到。
# 计算每个原子的logP和MR值：rdMolDescriptors._CalcCrippenContribs
# 返回结果是每个原子logP和MR元组的列表

from rdkit.Chem import rdMolDescriptors
contribs = rdMolDescriptors._CalcCrippenContribs(mol)
contribs[:3]


fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, [x for x,y in contribs])
fig.savefig('figures/descriptor_visulize3.png', bbox_inches='tight')





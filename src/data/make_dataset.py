# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import click
import logging
from pathlib import Path
import numpy as np
import joblib as jl 
import pandas as pd
import rdkit 
from rdkit.Chem import (MACCSkeys, RDKFingerprint, Descriptors, AllChem, rdchem, rdmolops)
from MolAnalysis.descriptors import CustomDescriptors

def checkIfSulphur(mol): 
    """ Checks if  an rdkit molecule has sulphur atoms.
    """
    functional_group = rdkit.Chem.MolFromSmiles("S")
    matches = mol.GetSubstructMatches(functional_group)
    return len(matches)

def read_xyz(path):

    atoms = []
    coordinates = []
    with open(path, 'r') as file:
        lines = file.readlines()
        n_atoms = int(lines[0])  
        smile = lines[n_atoms + 3].split()[0]        
    return smile

def process_wang(input_filepath, output_filepath):
    """ "Reads wang 2017 SMILES and saves them in a joblib dump"""
    wang = pd.read_csv(input_filepath+"wang/Supporting information CSV.csv")
    columns=wang.iloc[0]
    wang.rename(columns=wang.iloc[0], inplace=True)
    jl.dump(wang.iloc[1: , :], output_filepath+"wang_data.dump", compress=0, protocol=None, cache_size=None)

def process_gecko(input_filepath, output_filepath):
    geckofull = pd.read_csv(input_filepath+'gecko/Gecko_All_SMILES_full.txt')
    geckofull.drop('Unnamed: 0',axis=1,inplace=True)
    jl.dump(geckofull,output_filepath+'gecko_full.dump')

def process_quinones(input_filepath, output_filepath):
    quinones = pd.read_csv(input_filepath+'quinones/tabor_set_for_publ.csv')
    quinones.rename(columns={'smiles':'SMILES'}, inplace=True)
    CustomDescriptors(quinones)
    isS = pd.DataFrame(quinones[quinones.columns[quinones.columns.get_level_values(1)== 'mol']].values)[0].apply(lambda x: checkIfSulphur(x))
    jl.dump(pd.DataFrame(quinones.droplevel(0,axis=1)['SMILES'][isS==0]), output_filepath+'tabor_nosulf.dump')
    jl.dump(quinones.droplevel(0,axis=1), output_filepath+'tabor_all.dump')

def process_nabladft(input_filepath, output_filepath):
    """Reads presorted and unique nablaDFT data"""
    infile=input_filepath+'nablaDFT/summary.csv'
    result_smiles = subprocess.run(['/scratch/work/sandsth2/Projects/atm_datasets_manuscript/src/data/get_nabla_dft_smiles.sh', infile], stdout=subprocess.PIPE, text=True)
    pd.DataFrame(pd.DataFrame(result_smiles.stdout.splitlines()[1:])).to_csv(input_filepath+'nablaDFT/SMILES-sorted_unique.csv',header=None, index=False)
    nablaDFT = pd.read_csv(input_filepath+'nablaDFT/SMILES-sorted_unique.csv', header=None)
    nablaDFT.rename({0:'SMILES'},axis=1, inplace=True)
    jl.dump(nablaDFT, output_filepath+'nablaDFT.dump')

def process_qm9(input_filepath, output_filepath):
    """Read xyz files belonging to QM9 dataset and save their SMILES in a list"""
    XYZDIR = input_filepath+'qm9/xyz'
    data = []
    smiles = []
    properties = []
    IDs = []
    for file in os.listdir(XYZDIR):
        path = os.path.join(XYZDIR, file)
        smile = read_xyz(path)
        ID = int(file.split("_")[1].split(".")[0])
        smiles.append(smile) 
        IDs.append(ID)
    qm9_data = pd.DataFrame()
    qm9_data['ID'] = pd.DataFrame(IDs)
    qm9_data['SMILES'] = pd.DataFrame(smiles)
    jl.dump(qm9_data,output_filepath+"qm9_data.dump")

def process_all_mona(input_filepath, output_filepath):
    """Reads exp. mona data (preprocessed) and removes problematic smiles 
    that give error when coverting to rdkit molecule structure."""
    infile=input_filepath+'mona/all/MoNA-export-All_Spectra.msp'
    result_smiles = subprocess.run(['/scratch/work/sandsth2/Projects/atm_datasets_manuscript/src/data/get_smiles_from_mona_msp.sh', infile], stdout=subprocess.PIPE, text=True)
    pd.DataFrame(result_smiles.stdout.replace(',','').splitlines()).to_csv(output_filepath+"MoNA-export-All_Spectra.smiles.unique.-gt")
    MoNa_all = pd.read_csv(output_filepath+"MoNA-export-All_Spectra.smiles.unique.-gt",on_bad_lines='skip', low_memory=False);
    MoNa_all.rename(columns = {'0': "SMILES"}, inplace=True)
    MoNa_all = MoNa_all[MoNa_all['SMILES'].str.contains("InCh") == False]
    jl.dump(MoNa_all, output_filepath+"allmona_complete.dump")
    remove_smiles=['CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)',
    '(CC(=O)O3)O)C',
    'N/ACCC1(C(=O)NCNC1=O)c2ccccc2',
    'O=C1C2=C(C=C(C)C=C2O)OC3=CC(O)=CC(C(OC)=O)=C32',
    'OC1=CC=C(CC(C(NC(C(CC)C)C(OC(C(CCCCCCCCCC)C)CC(NC(C(NC(C(NC(C(NC2CCC(N)=O)=O)C)=O)C)=O)C(O)C)=O)=O)=O)NC2=O)C=C2',
    'OC1=CC(C(OC)=O)=C(OC2=CC(C)=CC(O)=C2C(O)=O)C(OC)=C2',
    'O=C([C@H](CC)C)O[C@H]1CCC=C2C1[C@@H](CC[C@@H](O)C[C@@H](O)CC(OC)=O)[C@@H](C)C=C3',
    'O=C(N(C(C=CC=C1)=C1C(N(C)[C@@]2([H])CC3=CC=CC=C3)=O)C2=N4)C5=C4C=CC=C6',
    'O=C(N[C@@H](CCCCCC(CC)=O)C(N[C@@H](CC1=CN(OC)C2=C1C=CC=C2)C3=O)=O)[C@@H]4N(C([C@H]([C@H](CC)C)N3)=O)CCCC5',
                  'CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)']
    MoNa_all = MoNa_all[MoNa_all['SMILES'].isin(remove_smiles) == False]
    bad_entries = []
    for j in MoNa_all.index: 
        try: 
            i = MoNa_all['SMILES'][j]
            mol = AllChem.MolFromSmiles(i)  
            fingerprint = MACCSkeys.GenMACCSKeys(mol)
            bit_fp = np.frombuffer(fingerprint.ToBitString().encode(), 'u1') - ord('0')
        except:
            bad_entries.append(j)
    jl.dump(MoNa_all.drop(bad_entries), output_filepath+'allmona_cleaned.dump')

def process_exp_mona(input_filepath, output_filepath):
    """Reads exp. mona data (preprocessed) and removes problematic smiles 
    that give error when coverting to rdkit molecule structure."""
    infile=input_filepath+'mona/exp/MoNA-export-Experimental_Spectra.msp'
    result_smiles = subprocess.run(['/scratch/work/sandsth2/Projects/atm_datasets_manuscript/src/data/get_smiles_from_mona_msp.sh', infile], stdout=subprocess.PIPE, text=True)
    pd.DataFrame(result_smiles.stdout.replace(',','').splitlines()).to_csv(output_filepath+"MoNA-export-Experimental_Spectra.smiles.unique.-gt")
    MoNa_exp = pd.read_csv(output_filepath+"MoNA-export-Experimental_Spectra.smiles.unique.-gt",on_bad_lines='skip', low_memory=False);
    MoNa_exp.rename(columns = {'0': "SMILES"}, inplace=True)
    MoNa_exp = MoNa_exp[MoNa_exp['SMILES'].str.contains("InCh") == False]
    jl.dump(MoNa_exp, output_filepath+"expmona_complete.dump")
    remove_smiles=['CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)',
    '(CC(=O)O3)O)C',
    'N/ACCC1(C(=O)NCNC1=O)c2ccccc2',
    'O=C1C2=C(C=C(C)C=C2O)OC3=CC(O)=CC(C(OC)=O)=C32',
    'OC1=CC=C(CC(C(NC(C(CC)C)C(OC(C(CCCCCCCCCC)C)CC(NC(C(NC(C(NC(C(NC2CCC(N)=O)=O)C)=O)C)=O)C(O)C)=O)=O)=O)NC2=O)C=C2',
    'OC1=CC(C(OC)=O)=C(OC2=CC(C)=CC(O)=C2C(O)=O)C(OC)=C2',
    'O=C([C@H](CC)C)O[C@H]1CCC=C2C1[C@@H](CC[C@@H](O)C[C@@H](O)CC(OC)=O)[C@@H](C)C=C3',
    'O=C(N(C(C=CC=C1)=C1C(N(C)[C@@]2([H])CC3=CC=CC=C3)=O)C2=N4)C5=C4C=CC=C6',
    'O=C(N[C@@H](CCCCCC(CC)=O)C(N[C@@H](CC1=CN(OC)C2=C1C=CC=C2)C3=O)=O)[C@@H]4N(C([C@H]([C@H](CC)C)N3)=O)CCCC5',
                  'CC(=O)CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)']
    MoNa_exp = MoNa_exp[MoNa_exp['SMILES'].isin(remove_smiles) == False]
    bad_entries = []
    for j in MoNa_exp.index: 
        try: 
            i = MoNa_exp['SMILES'][j]
            mol = AllChem.MolFromSmiles(i)  
            fingerprint = MACCSkeys.GenMACCSKeys(mol)
            bit_fp = np.frombuffer(fingerprint.ToBitString().encode(), 'u1') - ord('0')
        except:
            bad_entries.append(j)
    jl.dump(MoNa_exp.drop(bad_entries), output_filepath+'expmona_cleaned.dump')
    
def process_massbankeu(input_filepath, output_filepath):
    """Reads massbank 23.11 data """
    infile=input_filepath+'mb_eu/MassBank_NIST.msp'
    result_smiles = subprocess.run(['/scratch/work/sandsth2/Projects/atm_datasets_manuscript/src/data/get_smiles_from_mb_eu_msp.sh', infile], stdout=subprocess.PIPE, text=True)
    pd.DataFrame(result_smiles.stdout.replace(',','').splitlines()).to_csv(output_filepath+"mb_eu.smiles.unique.-gt")
    massbankeu  = pd.read_csv(output_filepath+"mb_eu.smiles.unique.-gt", on_bad_lines='skip', low_memory=False)
    massbankeu.rename(columns = {'0': "SMILES"}, inplace=True)
    jl.dump(massbankeu, output_filepath+"mb_eu_complete.dump") 
    bad_entries = []
    for j in massbankeu.index: 
        try: 
            i = massbankeu['SMILES'][j]
            mol = AllChem.MolFromSmiles(i)  
            fingerprint = MACCSkeys.GenMACCSKeys(mol)
            bit_fp = np.frombuffer(fingerprint.ToBitString().encode(), 'u1') - ord('0')
        except:
            bad_entries.append(j)
    massbankeu.drop(massbankeu.index[massbankeu["SMILES"].str.contains("\*") == True], inplace=True)
    jl.dump(massbankeu.drop(bad_entries), output_filepath+'mb_eu_cleaned.dump')
    
    
@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
    
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """
    logger = logging.getLogger(__name__)
    logger.info('making final nablaDFT data set from raw data')
    process_nabladft(input_filepath, output_filepath)
    logger.info('making final Wang data set from raw data')
    process_wang(input_filepath, output_filepath)
    logger.info('making final Gecko data set from raw data')
    process_gecko(input_filepath, output_filepath)
    logger.info('making final Quinones data set from raw data')
    process_quinones(input_filepath, output_filepath)
    logger.info('making final QM9 data set from raw data')
    process_qm9(input_filepath, output_filepath) 
    logger.info('making final MoNA data set from raw data')
    process_all_mona(input_filepath, output_filepath)
    logger.info('making final Massbank, Eu data set from raw data')
    process_massbankeu(input_filepath, output_filepath)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]
    
    main()

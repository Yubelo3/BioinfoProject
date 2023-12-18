from Bio.PDB import *
from Bio import SeqIO
import os
from tqdm import tqdm
import numpy as np
import warnings
from torch.utils.data import Dataset, DataLoader
import sys
sys.path.append("/home/Yubel/bioinfo/project/")
from data.tmalign import tm_score
import torch
from torch.nn.utils.rnn import pad_sequence


def pdbfile2fasta(path):
    with open(path,"r") as f:
        for record in SeqIO.parse(f,"pdb-atom"):
            return str(record.seq)

def fasta2int(fasta):
    ret=[]
    for c in fasta:
        assert 'A'<=c<='Z'
        ret.append(ord(c)-ord('A'))
    return ret

def load_database(target_file="data/database_full.npz",db_folder="pdb/"):
    '''
    如果target_file已经存在, 直接从target_file里读
    如果target_file不存在, 从db_folder中parse出需要的信息放在target_file里并返回
    '''
    if os.path.exists(target_file):
        res=np.load(target_file,allow_pickle=True)
        return res["seq"],res["endpoints"],res["filepaths"]
    
    seq,endpoints,filepaths=[],[0],[]  # 待存取信息
    # bar=tqdm(os.listdir(db_folder))
    bar=os.listdir(db_folder)
    process_cnt=0
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        for folder in bar:
            for file in os.listdir(db_folder+folder):
                filepath=folder+"/"+file
                try:
                    embedded_protein=fasta2int(pdbfile2fasta(db_folder+filepath))
                    endpoints.append(endpoints[-1]+len(embedded_protein))
                    seq.extend(embedded_protein)
                    filepaths.append(filepath)
                except:
                    print("invalid: "+filepath)
                process_cnt+=1
                # bar.set_description(f"{process_cnt}/107160 || {process_cnt/107160*100:.4f}%")
            # if process_cnt>=300:
            #     break
    np.savez_compressed(target_file,seq=seq,endpoints=endpoints,filepaths=filepaths)
    return seq,endpoints,filepaths

class PDBDataset(Dataset):
    def __init__(self,seq,endpoints,filepaths) -> None:
        super().__init__()
        self.seq=seq
        self.endpoints=endpoints
        self.filepaths=filepaths
    
    def __getitem__(self, index):
        return (torch.LongTensor(self.seq[endpoints[index]:endpoints[index+1]]),
                self.filepaths[index])

    def __len__(self):
        return len(self.filepaths)

def collate_fn(batch):
    n=len(batch)
    tm_scores=torch.zeros(n,n)
    for i in range(n):
        for j in range(n):
            tm_scores[i,j]=tm_score("pdb/"+batch[i][1],"pdb/"+batch[j][1])
    seq=pad_sequence([x[0] for x in batch],batch_first=True,padding_value=ord('Z')-ord('A')+1)
    return seq, tm_scores


if __name__=="__main__":

    seq,endpoints,filepaths=load_database("data/database_full.npz")
    # print(seq[:10])
    # print(endpoints[0])
    # print(endpoints[1])
    # print(endpoints[2])
    # print(filepaths[0])
    # dataset=PDBDataset(seq,endpoints,filepaths)
    # loader=DataLoader(dataset,batch_size=16,shuffle=True,collate_fn=collate_fn)
    # for x in loader:
    #     print(x[0].shape)
    #     print(x[1].shape)




# for folder in bar:
#     for file in os.listdir(prefix+folder):
        
#         parser = PDBParser(PERMISSIVE = True, QUIET = True) 
#         data = parser.get_structure("whatever","pdb/ZO/1ZO0-A.pdb")
#         print(data.seq)

        # dataset.append(1)
        # bar.set_description(str(len(dataset)))
        # print(data.header.keys())
        # models=data.get_models()
        # for model in models:
        #     chains=model.get_chains()  # 构成蛋白质的链
        #     for chain in chains:
        #         residues=chain.get_residues()  # 构成链的氨基酸残基
        #         for residue in residues:
        #             atoms=residue.get_atoms()  # 构成残基的原子
        #             for atom in atoms:
        #                 print(atom)
        #                 print(atom.get_vector())  # 原子的三维坐标
        # break


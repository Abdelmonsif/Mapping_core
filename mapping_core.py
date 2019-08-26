import pandas as pd
import argparse
import time

def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-inputcsv',default= '/home/monsif/Downloads/Mapping_core-master/pdb_aligning_results',required=True,help='aligned SNPs input file')
    parser.add_argument('-output',default= '/home/monsif/Downloads/Mapping_core-master/final_aligning',required=True,help='output of filtering')
    parser.add_argument('-highest',default= '/home/monsif/Downloads/Mapping_core-master/highest_pdb',required=True,help='highest pdb after filteration')
    parser.add_argument('-core',default= '/home/monsif/Downloads/Mapping_core-master/table_core',required=True,help='core table')
    return parser.parse_args()




def filtering(inputcsv,output,highest):
    df = pd.read_csv(inputcsv, delim_whitespace=True, names=['protein', 'SNP_ID', 'wt_codon', 'mu_codon', 'Mutant', 'Disorder', 'Condifdence', 'pdb', 'Position', 'Wild_type', 'PDB_position' , 'ID'])
    df1 = df[df.ID != '-']
    df1.to_csv(output, index=False, sep='\t')
    df2 = df1.groupby(['protein', 'SNP_ID'], sort=False, as_index=False).first()
    df2.to_csv(highest, index=False, sep='\t')




def mapping(inputcsv,core):
    df2 = pd.read_csv(highest, delim_whitespace=True)
    print('reading the file')
    with open(core, 'w') as h:
        for protein, SNP, position, wt_codon, mu_codon, wt, mutant, pdb, disorder, confidence, ID, pdb_position in zip(df2.protein, df2.SNP_ID, df2.Position, df2.wt_codon, df2.mu_codon,df2.Wild_type, df2.Mutant, df2.pdb, df2.Disorder, df2.Condifdence, df2.ID, df2.PDB_position):
            with open('/home/monsif/Downloads/Mapping_core-master/core/%s' %pdb, 'r') as f:
                lol = [x.replace('\n','') for x in f.readlines()]
                if str(position) in lol:
                    h.write(protein + '\t' + SNP + '\t' + str(position) + '\t' + wt + '\t' + mutant + '\t' + wt_codon + '\t' + mu_codon + '\t' + disorder + '\t' + str(confidence) + '\t' + pdb + '\t' + 'Not_IMP'+ '\t' + str(pdb_position) + '\t' + ID + '\t'  + '1' + '\t' + '0' + '\t' + '0' + '\n')
                else:
                    h.write(protein + '\t' + SNP + '\t' + str(position) + '\t' + wt + '\t' + mutant + '\t' + wt_codon + '\t' + mu_codon + '\t' + disorder + '\t' + str(confidence) + '\t' + pdb + '\t' + 'Not_IMP'+ '\t' + str(pdb_position) + '\t' + ID + '\t'  + '0' + '\t' + '0' + '\t' + '0' + '\n')
    h.close()




if __name__ == "__main__":
    args = getArgs()
    filtering = filtering(args.inputcsv,args.output,args.highest,args.core)
    #input = 
    
    
    start = time.time()
    end = time.time()
    print ('time elapsed:' + str(end - start))

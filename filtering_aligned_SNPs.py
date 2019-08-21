import pandas as pd

def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-input',default= '/home/monsif/Downloads/Mapping_core-master/pdb_aligning_results',required=True,help='aligned SNPs input file')
    parser.add_argument('-output',default= '/home/monsif/Downloads/Mapping_core-master/final_aligning',required=True,help='output of filtering')
    parser.add_argument('-highest',default= '/home/monsif/Downloads/Mapping_core-master/highest_pdb',required=True,help='highest pdb after filteration')
    return parser.parse_args()


def filter(input,output,highest):
    df = pd.read_csv(input, delim_whitespace=True, names=['protein', 'SNP_ID', 'wt_codon', 'mu_codon', 'Wild_type', 'Mutant', 'Position', 'Disorder', 'Condifdence', 'pdb', 'PDB_position' , 'ID'])
    df1 = df[df.ID != '-']
    df1.to_csv(output, index=False, sep='\t')
    df2 = df1.groupby(['protein', 'SNP_ID'], sort=False, as_index=False).first()
    df2.to_csv(highest, index=False, sep='\t')

if __name__ == "__main__":
    args = getArgs()
    filter = filter(args.input,args.output,args.highest)
    #input = 
    
    
    start = time.time()
    end = time.time()
    print ('time elapsed:' + str(end - start))

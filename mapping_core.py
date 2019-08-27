import pandas as pd
import argparse
import time

def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-inputcsv',default= '/home/monsif/Downloads/Mapping_core-master/pdb_aligning_results',required=True,help='aligned SNPs input file')
    parser.add_argument('-output',default= '/home/monsif/Downloads/Mapping_core-master/final_aligning',required=True,help='output of filtering')
    parser.add_argument('-highest',default= '/home/monsif/Downloads/Mapping_core-master/highest_pdb',required=True,help='highest pdb after filteration')
    parser.add_argument('-tablecore',default= '/home/monsif/Downloads/Mapping_core-master/table_core',required=True,help='core table')
    parser.add_argument('-core',default= '/home/monsif/Downloads/Mapping_core-master/core',required=True,help='core SNPs')
    parser.add_argument('-surface',default= '/home/monsif/Downloads/Mapping_core-master/surface',required=True,help='surface SNPs')
    parser.add_argument('-tableinterface',default= '/home/monsif/Downloads/Mapping_core-master/table_for_interface',required=True,help='interface table')
    return parser.parse_args()




def filtering(inputcsv,output,highest):
    df = pd.read_csv(inputcsv, delim_whitespace=True, names=['protein', 'SNP_ID', 'wt_codon', 'mu_codon', 'mutant', 'disorder', 'confidence', 'pdb', 'position', 'wildtype', 'pdb_position' , 'pdb_aa'])
    df1 = df[df.ID != '-']
    df1.to_csv(output, index=False, sep='\t')
    df2 = df1.groupby(['protein', 'SNP_ID'], sort=False, as_index=False).first()
    df2.to_csv(highest, index=False, sep='\t')




def mapping(highest,tablecore):
    df2 = pd.read_csv(highest, delim_whitespace=True)
    with open(tablecore, 'w') as h:
        for protein, SNP, position, wt_codon, mu_codon, wt, mutant, pdb, disorder, confidence, ID, pdb_position in zip(df2.protein, df2.SNP_ID, df2.Position, df2.wt_codon, df2.mu_codon, df2.Wild_type, df2.Mutant, df2.pdb, df2.disorder, df2.Confidence, df2.ID, df2.PDB_position):
            with open('/home/monsif/Downloads/Mapping_core-master/core/%s' %pdb, 'r') as f:
                core_res = [x.replace('\n','') for x in f.readlines()]
                if str(position) in core_res:
                    h.write(protein + '\t' + SNP + '\t' + str(position) + '\t' + wt + '\t' + mutant + '\t' + wt_codon + '\t' + mu_codon + '\t' + disorder + '\t' + str(confidence) + '\t' + pdb + '\t' + 'Not_IMP'+ '\t' + str(pdb_position) + '\t' + ID + '\t'  + '1' + '\t' + '0' + '\t' + '0' + '\n')
                else:
                    h.write(protein + '\t' + SNP + '\t' + str(position) + '\t' + wt + '\t' + mutant + '\t' + wt_codon + '\t' + mu_codon + '\t' + disorder + '\t' + str(confidence) + '\t' + pdb + '\t' + 'Not_IMP'+ '\t' + str(pdb_position) + '\t' + ID + '\t'  + '0' + '\t' + '0' + '\t' + '0' + '\n')
    h.close()
    
    


def interface_table(output,tablecore,tableinterface,core,surface):
    df3 = pd.read_csv(tablecore, delim_whitespace=True)
    df4 = df3[df3.core == 1]
    df4.to_csv(core, index=False, sep='\t')
    df5 = df3[df3.core == 1]
    df5.to_csv(surface, index=False, sep='\t')

    df = pd.read_csv(output, delim_whitespace=True)
    df1 = pd.read_csv(surface, delim_whitespace=True, usecols=['protein', 'SNP_ID', 'position', 'wildtype', 'mutant', 'disorder', 'confidence'])
    df2 = df1.merge(df, how = 'inner', on = ['protein', 'SNP_ID', 'position', 'wt_codon', 'mu_codon', 'disorder', 'confidence'])
    df2.to_csv(tableinterface, index=False, sep='\t')


def saving_interface_residues(tableinterface):
    df = pd.read_csv(tableinterface, delim_whitespace=True)
    chain_lst = df.pdb.str.get(4)
    for pdb,chain in zip(df.pdb, chain_lst):
        ppi_lst = [i.strip("\n").strip(" ") for i in list(open('/home/monsif/Downloads/Mapping_core-master/grep_interfaces/%s' %pdb))]
        if len(ppi_lst) == 0:
            pass
        else:
            for i in ppi_lst:
                df3 = pd.read_csv('/home/monsif/Downloads/Mapping_core-master/grep_interfaces/%s' %i, delim_whitespace=True, names=['atom','tw','wda', 'ldgme', 'chain', 'residue', 'sds', 'sdfw', 'sfww', 'fes', 'sdwef'])
                df4 = df3[df3['chain'] == chain]
                a7a = df4['residue']
                a7a.to_csv('/home/monsif/Downloads/Mapping_core-master/interface_residues/%s' %i, index=False, sep= '\t')



if __name__ == "__main__":
    args = getArgs()
    filtering = filtering(args.inputcsv,args.output,args.highest)
    mapping(args.highest,args.tablecore)
    interface_table(args.output,args.tablecore,args.tableinterface,args.core,args.surface)
    saving_interface_residues(args.tableinterface)
    #input = 
    
    
    start = time.time()
    end = time.time()
    print ('time elapsed:' + str(end - start))

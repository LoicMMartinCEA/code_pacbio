#!/usr/bin/env python3 
import pprint
import copy
from pickle import TRUE
import pandas as pd
import numpy as np

import argparse
import os,re,sys
import pickle

d3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
         'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
         'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
         'TRP':'W','TYR':'Y','TER':'Z'}
d1to3 = {}
for key,val in d3to1.items():
    d1to3[val] = key

restriction_enzymes = dict()
# restriction_enzymes['enzyme_name']  = (["site5'->3'","reverse5'->3'"], number of occurrences tolerated)
restriction_enzymes['BsaI']  = (['GGTCTC','GAGACC'], 0)
restriction_enzymes['NotI']  = (['GCGGCCGC'], 0) #palindromic
restriction_enzymes['Acc65I']= (['GGTACC'], 0) #palindromic
restriction_enzymes['KpnI']  = (['GGTACC'], 0) #palindromic
restriction_enzymes['AscI']  = (['GGCGCGCC'], 0) #palindromic
restriction_enzymes['EarI']  = (['CTCTTC','GAAGAG'], 0)
restriction_enzymes['BbsI']  = (['GAAGAC','GTCTTC'], 0)
restriction_enzymes['BsmBI'] = (['CGTCTC','GAGACG'], 0)
restriction_enzymes['NheI']  = (['GCTAGC'], 0) #palindromic
restriction_enzymes['SalI']  = (['GTCGAC'], 0) #palindromic
restriction_enzymes['PstI']  = (['CTGCAG'], 0) #palindromic
restriction_enzymes['SacI']  = (['GAGCTC'], 0) #palindromic
restriction_enzymes['SpeI']  = (['ACTAGT'], 0) #palindromic
restriction_enzymes['BamHI'] = (['GGATCC'], 0) #palindromic


dpert = {}
dpert['A'] = ['G','S','N','L','F','D','E','R','K']
dpert['C'] = ['A','G','S','L','F','D','E','R','K']
dpert['D'] = ['A','G','S','N','L','F','E','R','K']
dpert['E'] = ['A','G','S','Q','L','F','D','R','K']
dpert['F'] = ['A','G','S','N','L','M','E','R','K']
dpert['G'] = ['A','G','S','N','L','F','E','R','K']
dpert['H'] = ['A','G','S','L','F','D','E','R','K']
dpert['I'] = ['A','G','S','N','V','F','E','R','K']
dpert['K'] = ['A','G','S','N','L','F','D','E','R']
dpert['L'] = ['A','G','S','N','E','F','M','R','K']
dpert['M'] = ['A','G','S','N','L','F','E','R','K']
dpert['N'] = ['A','G','S','D','L','F','E','R','K']
dpert['P'] = ['A','G','S','N','L','F','E','R','K']
dpert['Q'] = ['A','G','S','N','L','F','E','R','K']
dpert['R'] = ['A','G','S','N','L','F','D','E','K']
dpert['S'] = ['A','G','T','N','L','F','E','R','K']
dpert['T'] = ['A','G','S','N','L','F','E','R','K']
dpert['V'] = ['A','G','S','N','L','F','E','R','K']
dpert['W'] = ['A','G','S','N','L','F','E','R','K']
dpert['Y'] = ['A','G','S','N','L','F','E','R','K']

dpert2fill20mut = dict()
dpert2fill20mut['A'] = ['V','Q','W','M']
dpert2fill20mut['C'] = ['V','Q','W','M']
dpert2fill20mut['D'] = ['V','Q','W','M']
dpert2fill20mut['E'] = ['V','Q','W','M']
dpert2fill20mut['F'] = ['D','Q','V','W']
dpert2fill20mut['G'] = ['D','Q','W','M','V']
dpert2fill20mut['H'] = ['V','Q','W','M']
dpert2fill20mut['I'] = ['D','Q','W','M']
dpert2fill20mut['K'] = ['V','Q','W','M']
dpert2fill20mut['L'] = ['D','Q','W','V']
dpert2fill20mut['M'] = ['D','Q','W','V']
dpert2fill20mut['N'] = ['V','Q','W','M']
dpert2fill20mut['P'] = ['D','Q','W','M']
dpert2fill20mut['Q'] = ['D','M','W','T','V']
dpert2fill20mut['R'] = ['V','Q','W','M']
dpert2fill20mut['S'] = ['D','Q','W','M']
dpert2fill20mut['T'] = ['D','Q','W','M']
dpert2fill20mut['V'] = ['D','Q','W','M']
dpert2fill20mut['W'] = ['D','Q','W','M']
dpert2fill20mut['Y'] = ['D','Q','W','M']


def read_fasta(file):
    str_seq = ""
    with open(file) as f:
        for ll in f.readlines():
            if ll[0] == ">":
                continue
            str_seq += ll.strip()
    return str_seq

def calc_deltanuc(old_codon,new_codon):
    diffnuc = 0
    for ii,nuc in enumerate(old_codon):
        if nuc != new_codon[ii]:
            diffnuc += 1
    return diffnuc

def nuc2aa(seq, dn2a, first_coding_nuc=1):
    aa_seq = ""
    for ii in range(0,len(seq[first_coding_nuc - 1:]),3):
        ii += first_coding_nuc - 1
        try:
            if d3to1[dn2a[seq[ii:ii+3]]] =='Z':
                break
            else:
                aa_seq += d3to1[dn2a[seq[ii:ii+3]]]
        except KeyError:
            print("Problem upon translation of your DNA coding sequence. "
                  "Either number of nucleotides is not a multiple of 3"
                  "or a codon cannot be translated")
            return None
    return aa_seq

def create_mutation_set(mut_file):
    set_mutres = set()
    p1 = re.compile("(\d+)")
    p2 = re.compile("(\d+)-(\d+)")
    with open(mut_file) as fmut:
        for ll in fmut.readlines():
            m1 = p1.findall(ll)
            m2 = p2.findall(ll)
            for resi in m1:
                set_mutres.add(int(resi))
            if len(m2) > 0:
                for res1,res2 in m2:
                    for resi in range(int(res1)+1,int(res2)):
                        set_mutres.add(int(resi))
    return set_mutres

def read_access(faccess):
    acc = open(faccess).readlines()
    dacc = {}
    count_aa = 1
    for ll in acc:
        sp = ll.split()
        if sp[0] == 'ACCESS':
            aa = d3to1[sp[2]]
            access = ll[62:67]
            count_aa += 1
            #print aa, count_aa, access
            dacc["%s%d"%(aa,count_aa)] = eval(access)
    return dacc

def check_restriction_site(nuc_seq):
    is_not_cut = True
    for enzyme_sites in restriction_enzymes.keys():
        number_of_instance_accepted = restriction_enzymes[enzyme_sites][1]
        for site in restriction_enzymes[enzyme_sites][0]:
            p = re.compile(site)
            m = p.findall(nuc_seq)
            if len(m) > number_of_instance_accepted:
                return "{}_{}".format(enzyme_sites, site)
    return is_not_cut

def create_python_script(target_name, fpathin, mutation_set):
    pymolfile = os.path.join(fpathin, "{}_mut_select.pml".format(target_name))
    str_mut = "+".join([str(x) for x in mutation_set])
    with open(pymolfile,"w") as pm:
        pm.write("cmd.select('mutations_{}','/{}///{}')\n".format(target_name, target_name, str_mut))
        pm.write("cmd.color('magenta', 'mutations_{}')\n".format(target_name))
    pm.close()
    return

def sample_mutants(target_name, fpathout, DNAseq, mut_index_set, thresh_abundance=12., index_start_codon = 1):

    twistlib_table_format = ['AA_wt', 'NUC_wt','PROTEIN_index','GCT', 'GCA', 'GCC', 'GCG', 'TGC', 'TGT', 'GAT', 'GAC', 'GAG', 'GAA', 'TTC', 'TTT', 'GGT', 'GGA', 'GGC', 'GGG', 'CAC', 'CAT', 'ATC', 'ATT', 'ATA', 'AAG', 'AAA', 'CTG', 'CTC', 'CTT', 'TTG', 'TTA', 'CTA', 'ATG', 'AAC', 'AAT', 'CCT', 'CCA', 'CCG', 'CCC', 'CAG', 'CAA', 'AGA', 'CGT', 'AGG', 'CGA', 'CGC', 'CGG', 'AGC', 'TCT', 'TCC', 'AGT', 'TCA', 'TCG', 'ACC', 'ACA', 'ACT', 'ACG', 'GTG', 'GTT', 'GTC', 'GTA', 'TGG', 'TAC', 'TAT', 'TGA', 'TAA', 'TAG']
    twistlib_AAconversion =  ['-','-','-'] + [dn2a[x] for x in twistlib_table_format[3:]]
    df = pd.DataFrame(columns = twistlib_table_format)
    print(twistlib_AAconversion)
    df.loc[0] = twistlib_AAconversion

    dabsort = {}

    for aa in da2n:
        if aa not in dabsort:
            dabsort[aa] = []
        for cod in da2n[aa]:
            if cod == "AAA": #quand c'est lysine 
                # Prevent long streches of AAA
                continue
            dabsort[aa].append([eval(da2n[aa][cod][0]),cod,aa]) # dic["AA"] : ["abudance","codon","AA"]
        if aa in ['LYS','PHE']:
            # To favour AAG instead of AAA for LYS and disfavour TTT for PHE
            dabsort[aa].sort()
        else:
            dabsort[aa].sort(reverse=True)
        dabsort[aa] = [x for x in dabsort[aa] if x[0] > thresh_abundance]#[:2]

    pprint.pprint(dabsort)

    d_alternative_mut = {}
    d_alternative_mut["E"] = [[62.83, 'GAT', 'ASP'], [37.17, 'GAC', 'ASP']]
    d_alternative_mut["S"] = [[46.83, 'ACC', 'THR'], [27.81, 'ACG', 'THR']]
    d_alternative_mut["N"] = [[66.6, 'CAG', 'GLN'], [33.4, 'CAA', 'GLN']]
    d_alternative_mut["R"] = [[75.44, 'AAA', 'LYS'], [24.56, 'AAG', 'LYS']]
    d_alternative_mut["L"] = [[51.2, 'ATT', 'ILE'], [44.37, 'ATC', 'ILE']]
    d_alternative_mut["A"] = [[43.17, 'GGC', 'GLY'], [32.91, 'GGT', 'GLY']]

    WT_seq = []
    dics = {}
    dic_mut = {}
    dic_mut['order'] = []
    dic_mut2addfor20mut = {}
    dic_mut2addfor20mut['order'] = []

    # If index_start_codon > 1 we can use the DNA tail in 5' to check for the absence of restriction site
    bias_start_index = index_start_codon - 1

    for sample in [target_name]:
        dics[sample] = {}
        dics[sample]['index'] = {}
        count_resindex = 0
        for count,cod in enumerate(range(bias_start_index,len(DNAseq),3)):
            codon = DNAseq[cod:cod+3]
            if d3to1[dn2a[codon]] == "Z":
                # STOP codon was reached
                break
            WT_seq.append(codon)

            index = count + 1

            dic_mut[count] = {}
            dic_mut[count]['order'] = []
            dic_mut['order'].append(count)

            dic_mut2addfor20mut[count] = {}
            dic_mut2addfor20mut[count]['order'] = []
            dic_mut2addfor20mut['order'].append(count)

            #print dn2a[codon],index
            print("codon:",codon,d3to1[dn2a[codon]],index)

            STOP_codon_frequency = 20

            resindex = index #"%s%d"%(d3to1[dn2a[codon]],index)
            if resindex in mut_index_set:
                count_resindex += 1
                list_all_mutants = copy.copy(dpert[d3to1[dn2a[codon]]])
                if count_resindex%STOP_codon_frequency == 0:
                    list_all_mutants.append('Z')
                    print("Introducing a STOP codon at",resindex)
                list_all_mutants.append(d3to1[dn2a[codon]])
                print(list_all_mutants)
                for kk, mutation in enumerate(list_all_mutants):

                    FLAG_MUT_OK = False
                    """
                    # CONDITION 1
                    """
                    aa3 = d1to3[mutation]
                    for abun_cod3 in dabsort[aa3]:
                        # Count how many nucleotide differ => should be above 1
                        deltanuc = calc_deltanuc(codon, abun_cod3[1])

                        if (deltanuc > 1 and abun_cod3[0] > thresh_abundance) or (deltanuc > 1 and abun_cod3[1] =='TTG') :
                            # TTG is a 'low abundant' 15% mutant and the second one for LEU...
                            # That's why we accepted the mutation as long as it had a deltanuc > 1
                            if mutation not in dic_mut[count]['order']:
                                dic_mut[count]['order'].append(mutation)
                                dic_mut[count][mutation] = []
                            dic_tmp = {}

                            header = ">mutAA:{}{}{} | mutnuc:{}->{} | nbmut:{} | freq: {:.1f} ".format(
                                d3to1[dn2a[codon]], resindex, mutation, codon, abun_cod3[1], deltanuc,abun_cod3[0])
                            print("Condition1: {}".format(header))
                            dic_tmp['real_mutation'] = mutation
                            dic_tmp['header'] = header
                            dic_tmp['new_codon'] = abun_cod3[1]
                            if d3to1[dn2a[codon]] == mutation:
                                dic_tmp['is_WT'] = True
                            else:
                                dic_tmp['is_WT'] = False
                            dic_mut[count][mutation].append(copy.copy(dic_tmp))
                            #nb_mut += 1
                            FLAG_MUT_OK = True
                    """    
                    # CONDITION 2: If condition 1 is not respected, we try to mutate into a closely related 
                    mutation which respects the conditions that more than one nucleotide changes are incorporated
                    """
                    if d3to1[aa3] in d_alternative_mut:
                        list_alt_mut = d_alternative_mut[d3to1[aa3]]
                        for abun_cod3 in list_alt_mut:
                            deltanuc = calc_deltanuc(codon, abun_cod3[1])
                            if deltanuc > 1 and abun_cod3[0] > thresh_abundance:
                                if mutation not in dic_mut[count]['order']:
                                    dic_mut[count]['order'].append(mutation)
                                    dic_mut[count][mutation] = []
                                dic_tmp = {}

                                alt_mutation = d3to1[abun_cod3[2]]

                                header = ">mutAA:{}{}{} (instead of {})| mutnuc:{}->{} | nbmut:{} | freq: {:.1f} ".format(
                                    d3to1[dn2a[codon]], resindex, alt_mutation, mutation, codon, abun_cod3[1], deltanuc,abun_cod3[0])
                                print("Condition2: {}".format(header))
                                dic_tmp['real_mutation'] = alt_mutation
                                dic_tmp['header'] = header
                                dic_tmp['new_codon'] = abun_cod3[1]
                                if d3to1[dn2a[codon]] == mutation:
                                    dic_tmp['is_WT'] = True
                                else:
                                    dic_tmp['is_WT'] = False
                                dic_mut[count][mutation].append(copy.copy(dic_tmp))

                                FLAG_MUT_OK = True

                    """    
                    # CONDITION 3 : We remove the condition on the number of differences
                    """
                    for abun_cod3 in dabsort[aa3]:
                        if abun_cod3[0] > thresh_abundance:
                            deltanuc = calc_deltanuc(codon, abun_cod3[1])
                            if deltanuc == 1:
                                if mutation not in dic_mut[count]['order']:
                                    dic_mut[count]['order'].append(mutation)
                                    dic_mut[count][mutation] = []
                                dic_tmp = {}

                                header = ">mutAA:{}{}{} | mutnuc:{}->{} | nbmut:{} | freq: {:.1f} ".format(
                                    d3to1[dn2a[codon]], resindex, mutation, codon, abun_cod3[1], deltanuc, abun_cod3[0])
                                print("Condition3: {}".format(header))
                                dic_tmp['real_mutation'] = mutation
                                dic_tmp['header'] = header
                                dic_tmp['new_codon'] = abun_cod3[1]
                                if d3to1[dn2a[codon]] == mutation:
                                    dic_tmp['is_WT'] = True
                                else:
                                    dic_tmp['is_WT'] = False
                                dic_mut[count][mutation].append(copy.copy(dic_tmp))

                                FLAG_MUT_OK = True

                    if not FLAG_MUT_OK:
                        print("################ Problem ###############",d3to1[dn2a[codon]], resindex,mutation)

                list_mutant2add20mut = copy.copy(dpert2fill20mut[d3to1[dn2a[codon]]])
                for kk, mutation in enumerate(list_mutant2add20mut):
                    aa3 = d1to3[mutation]
                    for abun_cod3 in dabsort[aa3]:
                        if abun_cod3[0] > thresh_abundance:
                            deltanuc = calc_deltanuc(codon, abun_cod3[1])
                            if deltanuc == 1:
                                continue
                            if mutation not in dic_mut2addfor20mut[count]['order']:
                                dic_mut2addfor20mut[count]['order'].append(mutation)
                                dic_mut2addfor20mut[count][mutation] = []
                            dic_tmp = {}
                            header = ">mutAA:{}{}{} | mutnuc:{}->{} | nbmut:{} | freq: {:.1f} ".format(
                                d3to1[dn2a[codon]], resindex, mutation, codon, abun_cod3[1], deltanuc, abun_cod3[0])
                            print("Additionals mutation to fill up to 20 mutants: {}".format(header))
                            dic_tmp['real_mutation'] = mutation
                            dic_tmp['header'] = header
                            dic_tmp['new_codon'] = abun_cod3[1]
                            dic_mut2addfor20mut[count][mutation].append(copy.copy(dic_tmp))

            dics[sample]['index'][index] = [codon,d3to1[dn2a[codon]],dn2a[codon]]

        dics[sample]['resi_last'] = index

    #pprint.pprint(dics)

    fout = open(os.path.join(fpathout,"{}_sequence_codons.fasta".format(target_name)),"w")

    ftwist = os.path.join(fpathout,"{}_all_mutated_sequences_for_twist.csv".format(target_name))

    max_total_mutations_per_position = 20

    # Now we have generated all the possible mutations in the right order.
    # We enumerate them, check for restriction sites (extending 5' for mutations in the first nucleotides).
    for pos in dic_mut['order']:
        df.loc[pos+1, 'AA_wt'] = dn2a[WT_seq[pos]]
        df.loc[pos+1, 'NUC_wt'] = WT_seq[pos]
        df.loc[pos+1, 'PROTEIN_index'] = pos+1
        set_mutated_codons = set()
        for mut in dic_mut[pos]['order']:
            max_number_of_solutions = 2
            count = 0
            count_accepted_solutions = 0
            for dposs_mut in dic_mut[pos][mut]:
                count += 1
                header = dposs_mut['header']
                new_codon = copy.copy(dposs_mut['new_codon'])
                if dposs_mut['is_WT']:
                    max_number_of_solutions = 1

                WT_seq_start = copy.copy(WT_seq[:pos])
                WT_seq_end = copy.copy(WT_seq[pos+1:])
                sequence_list =  WT_seq_start + [new_codon]  + WT_seq_end
                sequence = "".join(sequence_list) + "\n"

                is_not_cut = check_restriction_site(DNAseq[bias_start_index-6:bias_start_index]+sequence)

                if  is_not_cut is True:
                    count_accepted_solutions += 1
                    if count_accepted_solutions > max_number_of_solutions:
                        # Over the limit of synonymous mutations
                        break
                    df.loc[pos+1, new_codon] = new_codon
                    if new_codon not in set_mutated_codons:
                        fout.write(header)
                        fout.write(sequence)
                    set_mutated_codons.add(new_codon)
                else:
                    print("Found a case with an induced cleavage site {} which was corrected for ".format(is_not_cut),pos,mut)

        # Section to add the number of mutations required to have 20 mutations all in all
        if len(set_mutated_codons) < max_total_mutations_per_position:
            for mut in dic_mut2addfor20mut[pos]['order']:
                count_accepted_solutions = 0
                for dposs_mut in dic_mut2addfor20mut[pos][mut]:
                    count += 1
                    header = dposs_mut['header']
                    new_codon = copy.copy(dposs_mut['new_codon'])
                    WT_seq_start = copy.copy(WT_seq[:pos])
                    WT_seq_end = copy.copy(WT_seq[pos+1:])
                    sequence_list = WT_seq_start + [new_codon]  + WT_seq_end
                    sequence = "".join(sequence_list) + "\n"

                    is_not_cut = check_restriction_site(DNAseq[bias_start_index-6:bias_start_index] + sequence)

                    if is_not_cut is True:
                        count_accepted_solutions += 1
                        if count_accepted_solutions > max_number_of_solutions:
                            # Over the limit of synonymous mutations
                            break
                        print(">>> Added 1 mutation for position {}, ending up with {} mutations".format(pos, len(
                            set_mutated_codons)))
                        if len(set_mutated_codons) >= max_total_mutations_per_position:
                            # Over the number of mutations per positions
                            break

                        df.loc[pos+1, new_codon] = new_codon
                        if new_codon not in set_mutated_codons:
                            fout.write(header)
                            fout.write(sequence)
                        set_mutated_codons.add(new_codon)

                    else:
                        print("Found a case with an induced cleavage site {} which was corrected for ".format(is_not_cut),pos,mut)


    df = df.transpose()
    df.to_csv(ftwist)

    fout.close()


def load_pickles(fpath, list_pickle_files,verbose=False):
    list_dic = []
    for fpickle in list_pickle_files:
        print("Loading {}".format(fpickle))
        path_pickle = os.path.join(fpath, fpickle)
        try:
            dic_from_pickle = pickle.load(open(path_pickle))
        except TypeError:
            dic_from_pickle = pickle.load(open(path_pickle, 'rb'), encoding='latin1')
        if verbose:
            pprint.pprint(dic_from_pickle)
            print("############### end of {} ".format(fpickle))
        list_dic.append(dic_from_pickle)
    return list_dic


def cmdLine(debug=[None,None], verbose=False):
    USAGE = "A script to  blablbla.... \
        ./GenerateDMSLibrary.py -i machin.fasta -m truc.txt \
        If nothing is given it will use main code path" 
    #shit:USAGE = "A script to generate the table of mutants for synthesis by companies such as Twist, & co \n" \
    #shit:        "GenerateDMSMutantLibrary.py -i protein_a.pdb \n" #.... 
    

    parser = argparse.ArgumentParser(usage=USAGE)

    parser.add_argument('-i', '--seq_protein',
                        action='store',
                        dest='seq_protein',
                        help='Filename with the nucleotide sequence of the protein as fasta format',
                        default=debug[0]) 
    parser.add_argument('-m', '--mut_file',
                        action='store',
                        dest='mut_file',
                        help='Filename with the residue indexes ',
                        default=debug[1])
    return parser

if __name__=='__main__':
    # (name of the target_protein, first_coding_nuc in the AMPLICON sequence)
    #list_target = [("TNF",25),("VHH1",23),("VHH2",23),("VHH3",23)]#
    #list_target = [("TNF",25)]
    list_target = [("VHH2",23)]
    for target_protein,first_coding_nuc in list_target:

        if target_protein[:3] == "VHH":
            restriction_enzymes['PstI'] = (['CTGCAG'], 1)  # palindromic

        fpathin = "./"
        #fpathout = "./" doesn't seem to be defined anywhere... necessary ? 
        file_debug = os.path.join(fpathin,"{}_AMPLICON.fasta".format(target_protein))
        fmut_debug = os.path.join(fpathin,"{}_mut_index.txt".format(target_protein))
        # Location of the pickles
        fpath_pickles = "./"

        parser = cmdLine(debug=[file_debug, fmut_debug], verbose=False)
        options = parser.parse_args()

        pickle_files = ["nuc2aa.pickle", "aa2nuc.pickle", "nuc2abundance.pickle"]
        dn2a, da2n, dn2ab = load_pickles(fpath_pickles, pickle_files, verbose=False) #dn2ab is useless ? same as da2n ? 

        print("Reading the fasta nucleotide file : {}".format(options.seq_protein))

        nuc_seq = read_fasta(options.seq_protein)
        print("Nucleotide sequence: {}".format(nuc_seq))

        mutation_set = create_mutation_set(options.mut_file)
        print("Number of mutants: {}".format(len(mutation_set)))

        aa_seq = nuc2aa(nuc_seq,dn2a,first_coding_nuc)
        print("Translated Amino acid sequence in the Coding Frame: {}".format(aa_seq))

        sample_mutants(target_protein, fpathin, nuc_seq, mutation_set, index_start_codon=first_coding_nuc)

        create_python_script(target_protein,fpathin, mutation_set)
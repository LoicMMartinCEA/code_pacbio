#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

'''

### MODULES IMPORT ###
from os import listdir
from os.path import isfile, join
from datetime import datetime
import locale
import warnings
import utils
import pandas as pd
warnings.filterwarnings("ignore")

# American global settings
locale.setlocale(locale.LC_ALL, 'C')

### GLOBAL VARIABLES SECTION ###
OptFN = RefFastaFN = CCSFastaFN = OutFN = ""
nucmer = "nucmer"  # Path to nucmer executable - if in the system path then just the executable name is enough
show_snps = "show-snps" # Path to show-snps executable - if in the system path then just the executable name is enough
LogFN = 'treat_LRs.log'
UnExpMutFN = 'unexpected_mutations.txt'
Nthreads = 30
OptionsDic = {}
MMDic_glob = {}
CCSFastaDict = {}
RefsDic_glob = {}
ConstRefL = []

### LOGFILE AVAILABLE TO ANY FUNCTION ###
LogFile=open (LogFN, "w", encoding="utf-8")
UEMFile=open(UnExpMutFN, "w", encoding="utf-8")


### FUNCTIONs ###
# write CDSs in fasta format : to be used for the alignment phase
def write_CDSsFasta(RefsDic):
    for ref in RefsDic.keys():
        CDSsF=open(f"{ref}_CDS.fa", "w", encoding="utf-8")
        CDSsF.write(f">{ref}\n{RefsDic[ref][5]}\n")
        CDSsF.close()

# treat CCSs valitated by the function get_CCSDic
def take_it(ref, ccs_name, anchorL, anchorR, seq, DegRefBC, CCSsDic):  # take_it()
    CCSsDic[ccs_name][0] = ref # record reference
    BC = utils.find_between( seq, anchorL, anchorR ) # record BC
    if DegRefBC.match(BC):
        CCSsDic[ccs_name][1] = BC
        # write validated CCS to the corresponding fasta file
        CCSsF=open(f"{ref}_CCSs.fa", "a+", encoding="utf-8")
        CCSsF.write(f">{ccs_name}\n{seq}\n")
        CCSsF.close()
    return CCSsDic

# get CCS infos
def get_CCSDic(FastQFN,RefsDic):
    CCSsDic = {}
    FastQ = utils.get_file_content(FastQFN)
    ccs_name = None
    seq = ""
    #qual = "" 
    # Generate CCS fasta file for each reference and treat constant BCs
    n = 0
    for ref in RefsDic.keys(): 
        # this section is to ensure begining with an empty file that will appended in the next steps
        CCSsF=open(f"{ref}_CCSs.fa", "w", encoding="utf-8")
        CCSsF.close()
        BC = RefsDic[ref][6]
        Pattern = '[ACGT]+$'
        if (utils.re.compile(Pattern)).match(BC) : # if BC is not degenerated (constant sequence of the partner)
            print(f"Non-degenerated BC ({BC}) for reference {ref}. Assuming a constant gene sequence (no mutants) ! ")
            LogFile.write(f"Non-degenerated BC ({BC}) for reference {ref}. Assuming a constant gene sequence (no mutants) !\n")
            ccs_name = f"@artificial_CCS_{n}"
            while ccs_name in CCSsDic.keys():
                n += 1
                ccs_name = f"@artificial_CCS_{n}"
            CCSsDic[ccs_name] = [ref,BC,"WT",[],[]] #[ref, BC, Type, info_NT, info_AA]
            CCSsDic = take_it(ref, ccs_name, RefsDic[ref][1][0], RefsDic[ref][1][1], RefsDic[ref][4], \
                              utils.re.compile(utils.Nt2RE(RefsDic[ref][6])), CCSsDic)
            ConstRefL.append(ref)
        else:
            print(f"Degenerated BC ({BC}) for reference {ref}. Assuming that gene mutants are expected ! ")
            LogFile.write(f"Degenerated BC ({BC}) for reference {ref}. Assuming that gene mutants are expected !\n")

    # Analyze LR fastQ. Read fastq file once and for each read try to assign a corresponding reference using its index
    print(f"\n{datetime.now()} | Entering LR fastQ analysis phase")
    LogFile.write(f"\n{datetime.now()} | Entering LR fastQ analysis phase\n")
    i = 1
    for line in FastQ:
        if line[0] == "@" and line.split("/")[-1] == "ccs":
            ccs_name = line
            CCSsDic[ccs_name] = [None,None,"WT",[],[]] #[ref, BC, Type, info_NT, info_AA]
            i += 1
        else:
            if line != "+":
                if i == 2: 
                    seq = line # 2nd line of a fastq read = sequence
                    for ref in RefsDic.keys():
                        if ref not in ConstRefL:
                            if RefsDic[ref][7] in seq: # if a reference index is found in CCS sequence
                                CCSsDic = take_it(ref, ccs_name, RefsDic[ref][1][0], RefsDic[ref][1][1], seq,\
                                                   utils.re.compile(utils.Nt2RE(RefsDic[ref][6])), CCSsDic)
                                break
                            else:
                                seq = utils.rc_dna(seq)
                                if RefsDic[ref][7] in seq:
                                    CCSsDic = take_it(ref, ccs_name, RefsDic[ref][1][0], RefsDic[ref][1][1], seq,\
                                                       utils.re.compile(utils.Nt2RE(RefsDic[ref][6])), CCSsDic)
                                    break

            if i < 4:
                i += 1
            else:
                #qual = line # we can get the read quality if required in the future
                i = 1
    print(f"{datetime.now()} | Exiting LR fastQ analysis phase")
    LogFile.write(f"{datetime.now()} | Exiting LR fastQ analysis phase\n")

    print(f"\n{len(CCSsDic)} long reads will be analyzed now:")
    LogFile.write(f"\n{len(CCSsDic)} long reads will be analyzed now: \n")

    return CCSsDic # { ccs_name : [ref, BC, Type, [info_NT], [info_AA]]}

# get mutations reported in show-aligns output
def parse_SNPs(CCSsDic,SNPFN):
    try:
        SNPinfo = utils.get_file_content(SNPFN)[4:] # from the 5th line to the end
        a = 0
        PercLast = 0
        TotalRecs = len(SNPinfo)
        for line in SNPinfo:
            Fields = line.split("\t")
            QueryNuc = Fields[2]
            RefNuc = Fields[1]
            if QueryNuc in ['A', 'C', 'G', 'T', '.'] and RefNuc in ['A', 'C', 'G', 'T', '.']:
                RefPos = Fields[0]
                ccs_name = Fields[11]
                NTmuta = [RefNuc, RefPos, QueryNuc]
                CCSsDic[ccs_name][3].append(NTmuta)
                CCSsDic[ccs_name][2] = "Mut"
                a += 1
                PercNow = int((a * 100)/TotalRecs)
                if PercLast != PercNow:
                    print("  Gathering SNP data ... %i %s             \
                                                                " % (PercNow, "%"), end="\r", flush=True)
                    PercLast = PercNow
        print("  Gathering SNP data ... Done!      ")
        LogFile.write("  Gathering SNP data ... Done!      \n")
        return CCSsDic

    except ValueError:
        return "Could not get SNPs!"

def get_SNPs(CCSsDic,RefsDic):
    for ref in RefsDic.keys():
        Nucmer_prefix = f"{ref}_nucmer"
        DeltaFN = f"{Nucmer_prefix}.delta"
        SNPFN = f"{ref}_nucmer.snp"
        RefFaFN = f"{ref}_CDS.fa"
        CCSsFN = f"{ref}_CCSs.fa"
        # run Nucmer using CDS of each reference and the CCS that were validated for the same reference
        cmd = ("%s -t %i -p %s -l 6 -L 6 %s %s  2>/dev/null" % (nucmer, Nthreads,  Nucmer_prefix, RefFaFN, CCSsFN)) # Redirect errors to /dev/null; Nucmer version 4rc1
        print(f"Running Nucmer for {RefFaFN} (reference) vs {CCSsFN} (filtered CCSs) ...             ", end="\r", flush=True)
        utils.run(cmd, LogFile)
        print(f"Running Nucmer for {RefFaFN} (reference) vs {CCSsFN} (filtered CCSs) ... Done !       ")
        LogFile.write(f"Running Nucmer for {RefFaFN} (reference) vs {CCSsFN} (filtered CCSs) ... Done !       \n")
        # run show-snps and redirect the output to SNPfile
        cmd = ("%s %s %s 2>/dev/null" % (show_snps, "-q -T", DeltaFN)) # Redirect errors to /dev/null
        SNPFile=open(SNPFN, "w", encoding="utf-8")
        print("  Running show-snps ...           ", end="\r", flush=True)
        utils.run(cmd, SNPFile)
        SNPFile.close()
        print("  Running show-snps ... Done!      ")
        LogFile.write("  Running show-snps ... Done!      \n")
        # read SNPfile and put information into a dictionary and do some work
        CCSsDic = parse_SNPs(CCSsDic,SNPFN)

    return CCSsDic

# check if a mutated codon is expected
def check_expected_mutation(MMDic, ref, codonnbr, codon, ccs_name):
    Expected = False
    try:
        df = MMDic[ref]
        if df.loc[codonnbr][codon] == "PASS":
            Expected = True
    except KeyError:
        UEMFile.write(f"{ccs_name} ==> AA position {codonnbr} or/and codon '{codon}' does not correspond to an expected mutation of the reference '{ref}'\n")

    return Expected

# translate nt mutations and classify CCSs based on their types
def get_AAmut(CCSsDic,RefsDic,MMDic):
    LogFile.write("\nEntering mutation analysis phase ...\n")
    CDSsDic = {}

    for ref in RefsDic.keys():
        FastaDic = utils.get_fasta(f"{ref}_CDS.fa")
        CDSsDic[ref] = FastaDic[ref]

    for ccs_name in CCSsDic:
        if CCSsDic[ccs_name][2] != "WT":
            ref = CCSsDic[ccs_name][0]
            try:
                CDSseq = CDSsDic[ref]
            except KeyError:
                print("Some unexpected problem related to 'CDSseq = CDSsDic[ref]' occurred !")
            CodonDic = {}
            AAList = []
            for muta in CCSsDic[ccs_name][3]:
                Posnbr = int(muta[1])
                CDS_start = 1  # since all reads were converted to the same DNA strand and reference CDSs should be indicated from its 1st bases,
                #1 is the only possibility here !
                codon_info = [] # codon, codonnbr and positionnbr
                codonnbr = int(((Posnbr - CDS_start) / 3) + 1)
                CodonPosEval = ((Posnbr + 1) - CDS_start ) / 3
                if CodonPosEval %1 == 0: # base is in the 3rd position of the codon
                    codon_info = [CDSseq[Posnbr -3 : Posnbr], codonnbr, 3]
                if 0.4 < CodonPosEval %1 < 0.8: # base is in the 2nd position of the codon
                    codon_info = [CDSseq[Posnbr -2 : Posnbr + 1], codonnbr, 2]
                if 0 < CodonPosEval %1 < 0.4: # base is in the 1st position of the codon
                    codon_info = [CDSseq[Posnbr -1 : Posnbr + 2], codonnbr, 1]
                # generate new codon from mutation, translate and report the AA mutation.
                if codonnbr in CodonDic.keys():
                    new_codon = list(CodonDic[codonnbr][1])
                    new_codon[codon_info[2] - 1] = muta[-1]
                    new_codon = "".join(new_codon)
                    CodonDic[codonnbr][1] = new_codon
                else:
                    orig_codon = codon_info[0]
                    new_codon = list(orig_codon)
                    new_codon[codon_info[2] - 1] = muta[-1]
                    new_codon = "".join(new_codon)
                    CodonDic[codonnbr] = [orig_codon, new_codon]

            ExpectedList = []
            for codonnbr in CodonDic.keys():
                AAresult = [utils.translate_dna(CodonDic[codonnbr][0]), codonnbr, utils.translate_dna(CodonDic[codonnbr][1])]
                AAList.append(AAresult)
                Expected = check_expected_mutation(MMDic, ref, codonnbr, CodonDic[codonnbr][1], ccs_name)
                if Expected not in ExpectedList:
                    ExpectedList.append(Expected)

            CCSsDic[ccs_name][4] = AAList # record aa mutations
            if len(AAList) == 1:
                if len(ExpectedList) == 1 and ExpectedList[0] == True:
                    if AAList[0][0] == AAList[0][2]: # if wt synonymous mutation
                        CCSsDic[ccs_name][2] = 'SynOK'
                    else:
                        CCSsDic[ccs_name][2] = 'MutOK'

        else: # if WT
            if CCSsDic[ccs_name][0] == None: # no reference index found
                if CCSsDic[ccs_name][1] == None: # no valid BC found
                    CCSsDic[ccs_name][2] = "Unparsable"
                else: # valid BC found
                    CCSsDic[ccs_name][2] = "Unparsable_BCOK"
            else: # reference index found
                if CCSsDic[ccs_name][1] == None: # no valid BC found
                    CCSsDic[ccs_name][2] = "Unparsable_RefOK"
                else: # valid BC
                    CCSsDic[ccs_name][2] = "WtOK"
    
    UEMFile.close()		
    return CCSsDic

# BC analysis | CCSsDic structure = { ccs_name : [Ref, BC, Type, [NT_info], [AA_info]]}
def BCanal(CCSsDic, MinNbrCCSs, RatioDivCCSs, RefsDic):
    # gather BC infos by reference
    BCByRefDic = {}
    NewDic = {}
    for ccs_name in CCSsDic.keys():
        if CCSsDic[ccs_name][0] != None and CCSsDic[ccs_name][1] != None : 
            # treat only CCSs with an assigned reference and a valid BC
            Ref = CCSsDic[ccs_name][0]
            BC = CCSsDic[ccs_name][1]
            Type = CCSsDic[ccs_name][2]
            NT_info = ""
            AA_info = ""
            # compile nt and aa variants description
            if len(CCSsDic[ccs_name][3]) > 0: # treat NT mutations
                NT_mut_nbr = 1
                for muta in CCSsDic[ccs_name][3]: # iterate over NT mutations
                    if NT_mut_nbr == 1:
                        NT_info = f"{muta[0]}{muta[1]}{muta[2]}"
                    else:
                        NT_info = f"{NT_info},{muta[0]}{muta[1]}{muta[2]}"
                    NT_mut_nbr += 1
            if NT_info == "":
                NT_info = "WT"
            if len(CCSsDic[ccs_name][4]) > 0: # treat translation of NT mutations to AA
                AA_mut_nbr = 1
                for muta in CCSsDic[ccs_name][4]: # iterate over AA mutations
                    if AA_mut_nbr == 1:
                        AA_info = f"{muta[0]}{muta[1]}{muta[2]}"
                    else:
                        AA_info = f"{AA_info},{muta[0]}{muta[1]}{muta[2]}"
                    AA_mut_nbr += 1
            if AA_info == '':
                AA_info = 'WT'

            if Ref in BCByRefDic.keys(): # known ref
                if BC in BCByRefDic[Ref].keys(): # known BC
                    if NT_info in BCByRefDic[Ref][BC].keys(): # known NT variant
                        BCByRefDic[Ref][BC][NT_info][3] += 1 
                    else: # 1st observation of the NT_info
                        BCByRefDic[Ref][BC][NT_info] = [ AA_info, Type, ccs_name, 1 ]
                else: # 1st obsertation of the BC
                    BCByRefDic[Ref][BC] = { NT_info: [ AA_info, Type, ccs_name, 1 ] }
            else: # 1st observation of the reference
                BCByRefDic[Ref] = { BC : { NT_info: [ AA_info, Type, ccs_name, 1 ]}} # BCByRefDic Structure => BCByRefDic[Ref] = {BC : {NT_info: [AA_info, type, ccs_name, counter]}} 

    # get only reliable BCs to a new dictionary
    AllowedL = ['MutOK', 'SynOK', 'WtOK'] # list of interesting variants (expected)
    for Ref in BCByRefDic.keys():
        if Ref not in ConstRefL:
            for BC in BCByRefDic[Ref].keys():
                MasterL = []
                for NT_info in BCByRefDic[Ref][BC].keys():
                    MasterL.append([BCByRefDic[Ref][BC][NT_info][3], BCByRefDic[Ref][BC][NT_info][1], NT_info]) # [[variant_count, type, NT_info]]
                # verify if master / second >= RatioDivCCSs
                MasterL = sorted(MasterL, key=lambda x:x[0], reverse=True)
                if len(MasterL) > 1:
                    if MasterL[0][0] >= MinNbrCCSs: # filter by: min. number of CCSs supporting nt variant
                        Ratio = MasterL[0][0] / MasterL[1][0]
                        if Ratio >= RatioDivCCSs: # filter by ratio:  Master / second 
                            if MasterL[0][1] in  AllowedL: # filter by: type of variant
                                if Ref in NewDic.keys():
                                    NewDic[Ref][BC] = { MasterL[0][2] : BCByRefDic[Ref][BC][MasterL[0][2]] }
                                else:
                                    NewDic[Ref] = { BC : { MasterL[0][2] : BCByRefDic[Ref][BC][MasterL[0][2]] } }
                            else:  # master not variant of interest
                                continue
                        else: # failed ratio
                            continue
                    else: # failed min. CCSs
                        continue
                else: # Bc with only one variant
                    if MasterL[0][0] >= MinNbrCCSs : # filter by: min. number of CCSs supporting nt variant
                        if MasterL[0][1] in  AllowedL: # filter by: type of variant
                            if Ref in NewDic.keys():
                                NewDic[Ref][BC] = { MasterL[0][2] : BCByRefDic[Ref][BC][MasterL[0][2]] }
                            else:
                                NewDic[Ref] = { BC : { MasterL[0][2] : BCByRefDic[Ref][BC][MasterL[0][2]] } }
                        else:  # master not variant of interest
                            continue
                    else: # failed min. CCSs
                        continue

    # finally, arbitrarily include references with constant BC
    for Ref in ConstRefL:
        NewDic[Ref] = { RefsDic[Ref][6] : { "WT" : ['WT', 'WtOK', '@m54278_210612_182645/4194616/ccs', 1000] } }

    BCByRefDic = NewDic

    return BCByRefDic

### RUN SCRIPT ###
StartDT = datetime.now()
print(f"{str(StartDT)} ==> Started the LR run\n")
LogFile.write(f"{str(StartDT)} ==> Started the LR run\n\n")
## load required data
# the options
AllowedPars = ["RefSeqs", "CCSFastQ", "MutMatrixFolder", "Nthreads", "MinNbrCCSs", "RatioDivCCSs", "Clean"]
OptionsDic = utils.get_Options(OptionsDic, AllowedPars)
# the reference sequences
RefsDic_glob = utils.get_Refs(OptionsDic["RefSeqs"])
## actual pipeline
# write CDS references Fasta
write_CDSsFasta(RefsDic_glob)
# the ccs sequences (also write independent fasta files based on the corresponding reference - based on the index)
CCSsDic_glob = get_CCSDic(OptionsDic["CCSFastQ"],RefsDic_glob)
CCSsDic_glob = get_SNPs(CCSsDic_glob,RefsDic_glob)
# get the mutation matrices for libraries
TSVFiles = [f for f in [a for a in listdir(OptionsDic["MutMatrixFolder"]) if isfile(join(OptionsDic["MutMatrixFolder"], a))] \
            if (f.split('.')[-1] == 'tsv') and (f.split('.')[0] not in ConstRefL)] 
TSVFiles.sort()
for MMFN in TSVFiles:
    MMRef = "".join(MMFN.split(".")[:-1])
    MMDF = pd.read_csv(OptionsDic["MutMatrixFolder"]+MMFN, sep='\t')
    MMDic_glob[MMRef] = MMDF.set_index('AA_pos')
CCSsDic_glob = get_AAmut(CCSsDic_glob,RefsDic_glob,MMDic_glob)
# export CCSsDic as a pickle file for subsequent use
utils.data2pickleLR(CCSsDic_glob, "LR_CCSs_info.pkl") 
# export CCSsDic to TSV
CNames = ["reference", "barcode", "type", "nt_mutation", "aa_mutations"]
CCSsDF = pd.DataFrame.from_dict(CCSsDic_glob).transpose()
CCSsDF.index.name = "ccs"
CCSsDF.columns = CNames
CCSsDF = CCSsDF.sort_values(by=["reference", "barcode", "type"], ascending=False)
CCSsDF.to_csv("LR_CCSs_info.tsv", sep="\t")
del CCSsDF
# BC analysis
BCByRefDic_glob = BCanal(CCSsDic_glob, int(OptionsDic["MinNbrCCSs"]), \
                     float(OptionsDic["RatioDivCCSs"]), RefsDic_glob)
# export BCByRefDic by reference as a pickle file for subsequent use
for Ref_glob in BCByRefDic_glob.keys():
    utils.data2pickleLR(BCByRefDic_glob[Ref_glob], f"LR_BCs_{Ref_glob}.pkl")
# for each reference: export BCDic to TSV
print("\nBarcode analysis summary:")
LogFile.write("\nBarcode analysis summary:\n")
for Ref_glob in BCByRefDic_glob.keys():
    BCDF = pd.DataFrame.from_dict(BCByRefDic_glob[Ref_glob]).transpose()
    BCDF.index.name = "barcode"
    BCDF.to_csv(f"LR_BCs_{Ref_glob}.tsv", sep="\t")
    print(f"  {len(BCByRefDic_glob[Ref_glob])} barcodes were retained for reference {Ref_glob}")
    LogFile.write(f"  {len(BCByRefDic_glob[Ref_glob])} barcodes were retained for reference {Ref_glob}\n")

# clean intermediary files if specified in the options file
if OptionsDic["Clean"] == "yes" or OptionsDic["Clean"].lower() == "yes":
    utils.clean(list(RefsDic_glob.keys()), LogFile)
# final statements
EndDT = datetime.now()
print(f"\n{str(EndDT)} ==> Ended LR analysis")
LogFile.write(f"\n{str(EndDT)} ==> Ended LR analysis\n")
print(f"\nRun duration: {EndDT - StartDT}\n\n")
LogFile.write(f"Run duration: {EndDT - StartDT}")
LogFile.close()

exit(0)


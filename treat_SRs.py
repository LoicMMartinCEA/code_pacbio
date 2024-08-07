#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

'''

### MODULES IMPORT ###
from datetime import datetime
import locale
import sys
import warnings
import pickle
import statistics
import utils
import pandas as pd
warnings.filterwarnings("ignore")

# American global settings
locale.setlocale(locale.LC_ALL, 'C')

### GLOBAL VARIABLES SECTION
LogFN = 'treat_SRs.log'
OptionsDic_glob = {}
ZScoreTreshold = 1

### LOGFILE AVAILABLE TO ANY FUNCTION
try:
    LogFile=open(LogFN, "w", encoding='utf-8')
except FileNotFoundError:
    print(f"Could not open Log file: {LogFN}")
    sys.exit(1)

### FUNCTIONs ###
# get enrichment normalized by WT mean enrichment
def get_dBCDic_NormEnrich(dBCDic):  # dBCDic[Ccomplex] = {dBC : {B/A : [pA_NT_mut, pA_AA_mut, pA_type, pB_NT_mut, pB_ AA_mut, pB_type, count, freq] } }

    WT_types = ["WtOK", "SynOK"]

    for Ccomplex in dBCDic.keys():
        WTDic = {}
        WTAvEnrich = None
        WTAvSD = None
        print(f"\nWT analysis report for {Ccomplex} complex")
        LogFile.write(f"\nWT analysis report for {Ccomplex} complex\n")

        # first round is to gather WT related information
        for dBC in dBCDic[Ccomplex]:
            for BA in dBCDic[Ccomplex][dBC]:
                # filter WT vs WT (real WTOK of SynOK)
                if ( dBCDic[Ccomplex][dBC][BA][2] in WT_types ) and ( dBCDic[Ccomplex][dBC][BA][5] in WT_types ):
                    if dBC in WTDic.keys():
                        WTDic[dBC][BA] = dBCDic[Ccomplex][dBC][BA][7] # freq
                    else:
                        WTDic[dBC] = { BA : dBCDic[Ccomplex][dBC][BA][7] }# freq

         # second round is to calculate enrichment of WT dBCs
        WT_enrichL = []
        for dBC in WTDic.keys():
            if len(WTDic[dBC]) == 2: # if we observe the BC before and after selection
                enrich = WTDic[dBC]["After"] /  WTDic[dBC]["Before"] # frequency enrichment calculation
                WT_enrichL.append(enrich)

        # remove outliers
        print(f"  Found {len(WT_enrichL)} WT vs WT double barcodes in both populations (before/after)")
        LogFile.write(f"  Found {len(WT_enrichL)} WT vs WT double barcodes in both populations (before/after)\n")
        WT_enrichL = utils.remove_outiers(WT_enrichL, ZScoreTreshold) # list of values, ZScoreTreshold

        try:
            WTAvEnrich = statistics.mean(WT_enrichL)
            WTAvSD = statistics.stdev(WT_enrichL)
            print(f"  Kept {len(WT_enrichL)} WT vs WT double barcodes after outliers removal.\
                   WT average enrichment: {WTAvEnrich:.3} \u00B1 {WTAvSD:.3}")
            LogFile.write(f"  Kept {len(WT_enrichL)} WT vs WT double barcodes after outliers removal.\
                           WT average enrichment: {WTAvEnrich:.3} \u00B1 {WTAvSD:.3}")
        except ArithmeticError:
            print(f"  Could not calculate mean and sd from wild-type list:\n{WT_enrichL}")
            LogFile.write(f"  Could not calculate mean and sd from wild-type list:\n{WT_enrichL}\n")
        del WT_enrichL, WTDic

        # third round is calculate all Normalized scores (enrichments related to WT enrichment)
        NormEnrich = None
        for dBC in dBCDic[Ccomplex]:
            if len( dBCDic[Ccomplex][dBC] ) == 2: # if we observe the BC before and after selection
                enrich =  dBCDic[Ccomplex][dBC]["After"][7] /  dBCDic[Ccomplex][dBC]["Before"][7] # frequency enrichment calculation
                NormEnrich = enrich / WTAvEnrich
                dBCDic[Ccomplex][dBC]["NormEnrich"] = NormEnrich

    return dBCDic

# get dBC frequencies
def get_dBCDic_freqs(dBCDic, n_readsDic):
    for Ccomplex in dBCDic.keys():
        for dBC in dBCDic[Ccomplex]:
            for BA in dBCDic[Ccomplex][dBC]:
                dBCDic[Ccomplex][dBC][BA][7] = dBCDic[Ccomplex][dBC][BA][6] / n_readsDic[Ccomplex][BA]

    return dBCDic

# parse short read options
def get_dBCDic(OptionsDic, RefsDic):

    #retrieve information about the complexes to analyse and their corresponding FastQ files.
    #ComplexesL and BA_FastQL should have the same length in all dimensions
    ComplexesL = utils.double_split(OptionsDic["Partners"], ";", "~") 
    # Partners name and References name should be identical
    BA_FastQL = utils.double_split(OptionsDic["Before_AfterFQ"], ";", "~")
    dBCDic = {} # create double barcode dictionary
    n_readsDic = {} # dictionary to record the number of reads in samples (before/after)

    for i in range(len(ComplexesL)): # for each complex
        pA_name = ComplexesL[i][0]
        pA_anchors = RefsDic[pA_name][1]
        pA_LRDic = {}
        pB_name = ComplexesL[i][1]
        pB_anchors = RefsDic[pB_name][1]
        pB_LRDic = {}
        Ccomplex = f"{pA_name}_vs_{pB_name}" # current complex name
        print(f"Searching double barcodes for partners {Ccomplex}:")
        LogFile.write(f"Searching double barcodes for partners {Ccomplex}:\n")
        n_readsDic[Ccomplex] = {}
        with open(r"LR_BCs_{}.pkl".format(pA_name), "rb") as input_file:
            pA_LRDic = pickle.load(input_file)
        with open(r"LR_BCs_{}.pkl".format(pB_name), "rb") as input_file:
            pB_LRDic = pickle.load(input_file)
        for j in range(len(BA_FastQL[i])): # for each FastQ file
            FastQFN = BA_FastQL[i][j]
            FastQ = utils.get_file_content(FastQFN)
            seq = BA = ""
            if j == 0:
                BA = "Before"
            else:
                BA = "After"

            if BA not in n_readsDic[Ccomplex].keys():
                n_readsDic[Ccomplex][BA] = None
            n = 1
            n_reads = 0
            for line in FastQ:
                if line[0] == "@" and len(line.split(":")) == 7:
                    n += 1
                else:
                    if line != "+":
                        if n == 2:
                            seq = line # record sequence
                    if n < 4:
                        n += 1
                    else:
                        BC_pA = utils.find_between(seq, pA_anchors[0], pA_anchors[1])
                        BC_pB = utils.find_between(utils.rc_dna(seq), pB_anchors[0], pB_anchors[1])
                        if BC_pA in pA_LRDic.keys() and BC_pB in pB_LRDic.keys():
                            dBC = f"{BC_pA}:{BC_pB}"
                            pA_NT_mut = list(pA_LRDic[BC_pA].keys())[0]
                            pA_AA_mut = pA_LRDic[BC_pA][pA_NT_mut][0]
                            pA_type = pA_LRDic[BC_pA][pA_NT_mut][1]
                            pB_NT_mut = list(pB_LRDic[BC_pB].keys())[0]
                            pB_AA_mut = pB_LRDic[BC_pB][pB_NT_mut][0]
                            pB_type = pB_LRDic[BC_pB][pB_NT_mut][1]
                            if Ccomplex in dBCDic.keys():
                                if dBC in dBCDic[Ccomplex].keys():
                                    if BA in dBCDic[Ccomplex][dBC].keys():
                                        dBCDic[Ccomplex][dBC][BA][6] += 1
                                    else: # if 1st observation of the the sample (before or after)
                                        dBCDic[Ccomplex][dBC][BA] = \
                                            [ pA_NT_mut, pA_AA_mut, pA_type, pB_NT_mut, pB_AA_mut, pB_type, 1, None ]
                                else: # 1s observation of the dBC
                                    dBCDic[Ccomplex][dBC] = \
                                        { BA : [ pA_NT_mut, pA_AA_mut, pA_type, pB_NT_mut, pB_AA_mut, pB_type, 1, None ] }
                            else: # if 1st observation of the Ccomplex
                                dBCDic[Ccomplex] = \
                                    { dBC : { BA : [ pA_NT_mut, pA_AA_mut, pA_type, pB_NT_mut, pB_AA_mut, pB_type, 1, None ] } }
                        n = 1
                        n_reads += 1
            n_readsDic[Ccomplex][BA] = n_reads

        print(f"  {len(dBCDic[Ccomplex])} known and valid double barcodes were found in both samples (Before & After).")
        LogFile.write(\
            f"  {len(dBCDic[Ccomplex])} known and valid double barcodes were found in both samples (Before & After).\n")

    dBCDic = get_dBCDic_freqs(dBCDic, n_readsDic)
    dBCDic = get_dBCDic_NormEnrich(dBCDic)

    # dBCDic[Ccomplex] = {dBC : {B/A : [pA_NT_mut, pA_AA_mut, pA_type, pB_NT_mut, pB_ AA_mut, pB_type, count, freq] } }
    return dBCDic

# do an aa variant focused analysis
def get_AAVarStats(AAVarDic):
    for Ccomplex in AAVarDic.keys():
        for dMutAA in AAVarDic[Ccomplex].keys():
            AA_enrichL = []
            for IndBCData in AAVarDic[Ccomplex][dMutAA]['IndBCData']:
                AA_enrichL.append(IndBCData[1])
            AA_enrichL = utils.remove_outiers(AA_enrichL, ZScoreTreshold) # list of values, ZScoreTreshold
            if len(AA_enrichL) > 1:
                AAAvEnrich = statistics.mean(AA_enrichL)
                AAAvSD = statistics.stdev(AA_enrichL)
            else:
                AAAvEnrich = AA_enrichL[0]
                AAAvSD = None
            AAVarDic[Ccomplex][dMutAA]['Stats'] = [ AAAvEnrich, AAAvSD, len(AA_enrichL) ]

    return AAVarDic

# do an aa variant focused analysis
def get_AAVarDic(dBCDic):
    AAVarDic = {}
    WT_types = ["WT", "SynOK"]
    for Ccomplex in  dBCDic.keys():
        for dBC in dBCDic[Ccomplex].keys():
            if "NormEnrich" in dBCDic[Ccomplex][dBC].keys(): # if a normalized enrich could be calculated for the BC
                if Ccomplex in AAVarDic.keys():
                    mutA = dBCDic[Ccomplex][dBC]["Before"][1]
                    mutB = dBCDic[Ccomplex][dBC]["Before"][4]
                    if ( dBCDic[Ccomplex][dBC]["Before"][2] in WT_types ) :
                        mutA = "WT"
                    if ( dBCDic[Ccomplex][dBC]["Before"][5] in WT_types ) :
                        mutB = "WT"
                    dMutAA = f"{mutA}:{mutB}"
                    if dMutAA in AAVarDic[Ccomplex].keys():
                        AAVarDic[Ccomplex][dMutAA]['IndBCData'].append([dBC, dBCDic[Ccomplex][dBC]['NormEnrich']])
                    else:
                        AAVarDic[Ccomplex][dMutAA] = \
                            {'IndBCData': [[dBC, dBCDic[Ccomplex][dBC]['NormEnrich']]], 'Stats' : [None, None, None] }
                        # AAVarDic[Ccomplex][dMutAA] = {'IndBCData': [[dBC, NormEnrich]], 'Stats' : [Mean, SD, NbrBC]}
                else:
                    mutA = dBCDic[Ccomplex][dBC]["Before"][1]
                    mutB = dBCDic[Ccomplex][dBC]["Before"][4]
                    if ( dBCDic[Ccomplex][dBC]["Before"][2] in WT_types )  :
                        mutA = "WT"
                    if ( dBCDic[Ccomplex][dBC]["Before"][5] in WT_types ) :
                        mutB = "WT"
                    dMutAA = f"{mutA}:{mutB}"
                    AAVarDic[Ccomplex] = \
                          { dMutAA: {'IndBCData': [[dBC, dBCDic[Ccomplex][dBC]['NormEnrich']]], 'Stats' : [None, None, None] } }
                    # AAVarDic[Ccomplex]= {'IndBCData': [[dBC, NormEnrich]], 'Stats' (- outliers): [Mean enrichment, SD, Nbr of different dBCs ]} }

    del dBCDic
    AAVarDic = get_AAVarStats(AAVarDic)
    return AAVarDic

### RUN SCRIPT ###
StartDT = datetime.now()
print(f"{str(StartDT)} ==> Started the SR run\n")
LogFile.write(f"{str(StartDT)} ==> Started the SR run !\n\n")
## load required data
# the options
AllowedPars = ["RefSeqs", "Partners", "Before_AfterFQ", "Clean", "ZScoreTreshold"]
OptionsDic_glob = utils.get_Options(OptionsDic_glob, AllowedPars)
ZScoreTreshold = int(OptionsDic_glob["ZScoreTreshold"])
# the reference sequences
RefsDic_glob = utils.get_Refs(OptionsDic_glob["RefSeqs"])
# actully do the analysis
dBCDic_glob = get_dBCDic(OptionsDic_glob, RefsDic_glob) # get double barcode information
# for each commplex export dBCDic to .tsv and .pkl
for Ccomplex_glob in dBCDic_glob.keys():
    dBCDF = pd.DataFrame.from_dict(dBCDic_glob[Ccomplex_glob]).transpose()
    dBCDF.index.name = "dBC"
    dBCDF.to_csv(f"SR_dBC_{Ccomplex_glob}.tsv", sep="\t")
    utils.data2pickle(dBCDic_glob, Ccomplex_glob, f"SR_dBC_{Ccomplex_glob}.pkl")
AAVarDic_glob = get_AAVarDic(dBCDic_glob) # get AA variant focused analysis
# for each complex export AAVarDic_glob to .tsv and .pkl
for Ccomplex_glob in AAVarDic_glob.keys():
    AAvarDF = pd.DataFrame.from_dict(AAVarDic_glob[Ccomplex_glob]).transpose()
    AAvarDF.index.name = "AA_variant"
    AAvarDF.to_csv(f"SR_AAvars_{Ccomplex_glob}.tsv", sep="\t")
    utils.data2pickle(AAVarDic_glob, Ccomplex_glob, f"SR_AAvars_{Ccomplex_glob}.pkl")

# final statements
EndDT = datetime.now()
print(f"\n{str(EndDT)} ==> Ended SR analysis")
LogFile.write(f"\n\n{str(EndDT)} ==> Ended SR analysis\n")
print(f"\nRun duration: {EndDT - StartDT}")
LogFile.write(f"Run duration: {EndDT - StartDT}")

LogFile.close()

exit(0)


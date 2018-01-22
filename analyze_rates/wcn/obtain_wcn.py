"""
    Get WCN for all enzymes and see how inferred model affects relationship
    Note we need to only consider sites which are not gaps in the PDB reference of the alignment
"""
import os
from Bio import AlignIO

wcn_directory = "/Users/sjspielman/Research/enzyme_distance/enzyme_data/wcn_mono/"
data_directory =  "../data_preparation/enzymes/phy/"

my_enzymes = ["13PK_A", "1B93_A", "1CF2_O", "1E6E_A", "1G79_A", "1HV9_A", "1MQW_A", "1QAZ_A", "1STC_E", "1YON_A", "1A16_A", "1BF2_A", "1CTN_A", "1E7Q_A", "1G8F_A", "1IVH_A", "1NBA_A", "1QFM_A", "1TAH_B", "1YVE_I", "1A4L_A", "1BJO_A", "1CTT_A", "1EG7_A", "1GA8_A", "1J79_A", "1NIR_A", "1QI9_A", "1TLP_E", "1ZE1_A", "1A50_B", "1BQC_A", "1CZ1_A", "1EUU_A", "1GCB_A", "1JOF_A", "1NU3_A", "1QJE_A", "1TRK_B", "206L_A", "1A8H_A", "1BZC_A", "1DBT_A", "1EXP_A", "1GDH_A", "1K30_A", "1OAC_A", "1QRG_A", "1TYF_A", "2DOR_A", "1A8S_A", "1C0K_A", "1DEK_A", "1FF3_A", "1GP1_A", "1L1D_A", "1OH9_A", "1REQ_A", "1UAQ_A", "2JCW_A", "1AF7_A", "1C2T_A", "1DJ0_A", "1FPS_A", "1GPJ_A", "1L7N_A", "1OPM_A", "1S3I_A", "1UQT_A", "3CSM_A", "1AJ8_A", "1C82_A", "1DPG_A", "1FR8_A", "1GTP_A", "1LJL_A", "1OS7_A", "1S95_A", "1W0H_A", "3ECA_A", "1B73_A", "1CBG_A", "1E2A_A", "1G6T_A", "1H3I_A", "1MBB_A", "1PMA_B", "1SML_A", "1WGI_A", "3PVA_A", "1B8F_A", "1CEL_A", "1E3V_A", "1G72_A", "1HR6_B", "1MHY_D", "1POW_A", "1SMN_A", "1YCF_A", "7ATJ_A"]

outheader = "dataset,site,wcnSC\n"
outstring = ""
for dataset in my_enzymes:
   
    aln = AlignIO.read(data_directory + dataset + ".phy", "phylip-relaxed")
    for i in range(len(aln)):
        if "lcl|UniRef90" not in str(aln[i].id):
            ref = aln[i]
            break
    seq = str(ref.seq)
    indices = [] ## not gaps in reference PDB
    for i in range(len(seq)):
        if seq[i] == "-":
            continue
        else:
            indices.append(i)

    wcn_file = wcn_directory + dataset + "_wcn.csv"    
    with open(wcn_file, "r") as f:
        wcnraw = f.read().strip()
    wcnlines = wcnraw.split("\n")[1:]
    
    if( len(indices) != len(wcnlines)):
        continue
    for line in wcnlines[1:]:
        line2 = line.strip().split(",")
        
        site = str(indices[ int(line2[1])-1 ])
        wcnSC = line2[2]
        
        outstring += dataset + "," + site + "," + wcnSC + "\n"
    
outstring.strip()
with open("enzyme_wcnSC.csv", "w") as f:
    f.write(outheader + outstring)
    

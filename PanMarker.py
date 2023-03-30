import argparse
from functools import partial
import os
from multiprocessing.pool import Pool
from scipy.stats import levene, ttest_ind
from outliers import smirnov_grubbs as grubbs
from scipy.stats import pearsonr
import numpy as np

def get_sample(input_file):
    global sample_list
    sample_list = []
    samples = {}
    with open(input_file,"r") as f1:
        for line in f1:
            a = line.strip().split(".")[-1]
            sample = line.strip().replace("."+a,"")
            sample_list.append(sample)
            samples[sample] = {}
    genelist = []
    for sample in samples:
        with open(sample+"."+a,"r") as f2:
            for line in f2:
                if ">" in line:
                    geneid = line.strip().replace(">","")
                    genelist.append(geneid)
                    samples[sample][geneid] = ""
                    continue
                samples[sample][geneid] += line.strip()
    genelist = list(set(genelist))
    return samples,genelist

def write_genefasta(samples,genelist):
    folder = os.path.exists("genelist_file")
    if not folder:
        os.makedirs("genelist_file") 
    for geneid in genelist:
        with open("genelist_file/"+geneid.replace(">",""),"w") as f1:
            for id in samples:
                try:
                    f1.write(">"+id+"\n"+samples[id][geneid]+"\n")
                except:
                    print(id+" don't have "+geneid)

def aln(geneid) :
    in_aln = "aln_file/"+geneid+".aln"
    run = "mafft --adjustdirection genelist_file/"+geneid+" > aln_file/"+geneid+".aln"
    os.system(run)
    os.system("rm genelist_file/"+geneid)
    out_var = in_aln.replace("aln","var")
    try:
        extract_var(in_aln, out_var)
    except:
        os.system(run)
    try:
        extract_var(in_aln, out_var)
    except:
        print(in_aln+" do not have various")
    os.system("rm "+in_aln)  
    out_stat = out_var.replace("var","stat")
    stat_var_vs_exp(out_var,out_stat)
    os.system("rm "+out_var)   
    write(out_stat,out_file)        
            
####################

def load_aln(aln_db, in_aln):
    with open(in_aln, 'r') as fin:
        id = ""
        seq = ""
        for line in fin:
            if line[0] == '>':
                if seq != "":
                    aln_db[id] = seq
                id = line.strip().split()[0][1:].replace('_R_', '')
                seq = ""
            else:
                seq += line.strip().upper()
        aln_db[id] = seq


def get_consensus_seq(aln_db):
    consensus_seq = ""
    smp_list = sorted(aln_db.keys())
    for _ in range(len(aln_db[smp_list[0]])):
        cnt_db = {}
        for smp in smp_list:
            base = aln_db[smp][_]
            if base not in cnt_db:
                cnt_db[base] = 0
            cnt_db[base] += 1
        max_base = ""
        max_cnt = 0
        for base in cnt_db:
            if cnt_db[base] > max_cnt:
                max_cnt = cnt_db[base]
                max_base = base
        if max_cnt == len(smp_list):
            consensus_seq += max_base.upper()
        else:
            consensus_seq += max_base.lower()
    return consensus_seq


def classify_info(classified_info, full_info):
    last_info = []
    for _ in range(len(full_info)):
        cur_info = full_info[_]
        if cur_info[1] == '-':
            type = '<INS>'
        elif cur_info[2] == '-':
            type = '<DEL>'
        else:
            type = '<SNP>'
        if last_info == []:
            last_info = [type, cur_info[0], cur_info]
        else:
            if last_info[0] == '<SNP>':
                classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                last_info = [type, cur_info[0], cur_info]
            else:
                if type != last_info[0]:
                    classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                    last_info = [type, cur_info[0], cur_info]
                else:
                    if cur_info[0] != last_info[1]+1:
                        classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                        last_info = [type, cur_info[0], cur_info]
                    else:
                        if last_info[2][3] != cur_info[3]:
                            classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])
                            last_info = [type, cur_info[0], cur_info]
                        else:
                            last_info[1] = cur_info[0]
                            if cur_info[1] != '-':
                                last_info[2][1] += cur_info[1]
                            if cur_info[2] != '-':
                                last_info[2][2] += cur_info[2]
    classified_info.append([last_info[2][0], last_info[0], last_info[2][1], last_info[2][2], last_info[2][3]])





def extract_var(in_aln, out_stat):

    aln_db = {}
    load_aln(aln_db, in_aln)

    with open(out_stat, 'w') as fout:
        smps = []
        ref_seq = get_consensus_seq(aln_db)
        seq_len = len(ref_seq)
        for _ in aln_db:
            diff_cnt = 0
            for __ in range(seq_len):
                if aln_db[_][__] != ref_seq[__].upper():
                    diff_cnt += 1
            if diff_cnt < seq_len*.2:
                smps.append(_)
                continue
        #fout.write("#>Consensus\n#%s\n"%ref_seq)
        fout.write("#POS\tTYPE\tREF\tALT\t%s\n"%('\t'.join(smps)))
        if len(smps) < len(aln_db)*.5:
            return
        # Indetify each pos
        full_info = []
        for i in range(seq_len):
            info = []
            ref = ref_seq[i].upper()
            alt = ""
            alt_cnt = {}
            for smp in smps:
                base = aln_db[smp][i]
                if base != ref:
                    type = 1
                    if base not in alt_cnt:
                        alt_cnt[base] = 0
                    alt_cnt[base] += 1
                else:
                    type = 0
                info.append(type)
            max_cnt = 0
            for base in alt_cnt:
                if alt_cnt[base] > max_cnt:
                    alt = base
                    max_cnt = alt_cnt[base]
            if alt == "" or max_cnt <= 1:
                continue
            full_info.append([i, ref, alt, '\t'.join(map(str, info))])
        classified_info = []
        classify_info(classified_info, full_info)
        for info in classified_info:
            fout.write("%s\n"%('\t'.join(map(str, info))))
    


def stat_var_vs_exp(in_var, out_stat):
    in_type = args.type 
    if in_type == "cds" :
        in_exp = phe
        exp_db = {}
        with open(in_exp, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                exp_db[data[0]] = float(data[1]) #np.log2(float(data[1])+1)
    geneid = in_var.replace("var_file/","").replace(".var","")
    if in_type == "prm"  :
        prm = {}
        with open(exp,"r") as f:
            for line in f:
                if "Gene_id" in line:
                    geneids = line.strip().split()[1:]
                    for num in range(len(geneids)):
                        geneid = geneids[num]
                        prm[geneid] = {}
                else:
                    lines = line.strip().split()
                    sampleid = lines[0]
                    for num in range(len(geneids)):
                        geneid = geneids[num]
                        FPKM = float(lines[num+1])
                        prm[geneid][sampleid] = FPKM
        geneid = in_var.replace("var_file/","").replace(".var","")
        exp_db = prm[geneid]
    with open(in_var, 'r') as fin:
        with open(out_stat, 'w') as fout:
            idx_db = {}
            for line in fin:
                data = line.strip().split()
                if line[0] == '#':
                    fout.write("%s\tLevene_pvalue\tTtest_pvalue\tValid_ref_count\tValid_alt_count\n"%line.strip())
                    smp_list = []
                    for i in range(4, len(data)):
                        idx_db[i] = data[i]
                        if idx_db[i] not in exp_db:
                            continue
                        smp_list.append(idx_db[i])
                else:
                    ref_list = []
                    alt_list = []
                    for i in range(4, len(data)):
                        if idx_db[i] not in exp_db:
                            continue
                        if data[i] == '0':
                            ref_list.append(exp_db[idx_db[i]])
                        elif data[i] == '1':
                            alt_list.append(exp_db[idx_db[i]])
                    # remove outliers with grubbs test
                    if len(ref_list) < 2 or len(alt_list) < 2:
                        continue
                    ref_list = grubbs.test(np.array(ref_list), alpha=0.05)
                    alt_list = grubbs.test(np.array(alt_list), alpha=0.05)
                    if len(ref_list) == 1 or len(alt_list) == 1:
                        continue
                    levene_val = levene(ref_list, alt_list).pvalue
                    if levene_val > 0.05:
                        equal_var = True
                    else:
                        equal_var = False
                    #equal_var = True
                    t = ttest_ind(ref_list, alt_list, equal_var=equal_var)
                    t_pval = t.pvalue
                    #is_write = False
                    #if np.average(ref_list) > np.average(alt_list) and min(ref_list) > max(alt_list):
                    #    is_write = True
                    #elif np.average(ref_list) < np.average(alt_list) and max(ref_list) < min(alt_list):
                    #    is_write = True
                    if in_type == "cds":
                        if t_pval <= 0.05  and abs(person[geneid]) >= 0.3:#and is_write:
                            fout.write("%s\t%.4f\t%.30f\t%d\t%d\n"%(line.strip(), levene_val, t_pval, len(ref_list), len(alt_list)))
                    elif in_type == "prm":  
                        if t_pval <= 0.05 :
                            fout.write("%s\t%.4f\t%.30f\t%d\t%d\n"%(line.strip(), levene_val, t_pval, len(ref_list), len(alt_list)))
                    #fout.write("%s\t%.4f\t%.30f\t%d\t%d\n"%(line.strip(), levene_val, t_pval, len(ref_list), len(alt_list)))

######################
def trait(expfile,phefile):
    exp_dict = {}
    phe_dict = {}
    with open(expfile,"r") as f1:
        for line in f1:
            if "Gene_id" in line:
                geneids = line.strip().split()[1:]
                for num in range(len(geneids)):
                    geneid = geneids[num]
                    exp_dict[geneid] = {}
            else:
                lines = line.strip().split()
                sampleid = lines[0]
                for num in range(len(geneids)):
                    geneid = geneids[num]
                    FPKM = float(lines[num+1])
                    exp_dict[geneid][sampleid] = FPKM
    with open(phefile,"r") as f2:    
        for line in f2:
            lines = line.strip().split()
            sampleid = lines[0]
            phenum = lines[1]
            phe_dict[sampleid] = phenum
    Person_dict = {}
    for geneid in  exp_dict:
        exp_list = []
        phe_list = []
        for sampleid in exp_dict[geneid]:
            exp_list.append(float(exp_dict[geneid][sampleid]))
            phe_list.append(float(phe_dict[sampleid]))
        array1 = np.array(exp_list)
        array2 = np.array(phe_list)
        pvalue =np.corrcoef(array1,array2)[0][1]
        Person_dict[geneid] = pvalue
    return Person_dict


    
def write(out_stat,outfile):
    geneid = out_stat.replace("stat_file/","").replace(".stat","")
    with open(out_stat,"r") as f1:
        for line in f1:
            lines = line.strip().split()
            if "#" in line:
                list = lines[4:-4]
                continue
            with open(outfile+".result","a") as f2:
                f2.write(geneid+"\t")
                f2.write(lines[0]+"\t"+lines[1]+"\t"+lines[2]+"\t"+lines[3]+"\t")
                for id in sample_list:#list1.index
                    if id in list:
                        f2.write(lines[4+list.index(id)]+"\t")
                    else:
                        f2.write("NA\t")
                f2.write(lines[-4]+"\t"+lines[-3]+"\t"+lines[-2]+"\t"+lines[-1]+"\n")


######################


def extract_var_result(in_aln):
    out_var = in_aln.replace("aln","var")
    try:
        extract_var(in_aln, out_var)
    except:
        print(in_aln+" do not have variation")
    out_stat = out_var.replace("var","stat")
    try:
        stat_var_vs_exp(out_var,out_stat)
    except:
        print(out_var+" is filtered")
    os.system("rm "+out_var)


###########


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("-i","--inputfile",help="input cds files list",type=str,required=True)#需要参与比对的id文件 
parser.add_argument("-t","--threat",help="num of threats",type=int, default=10)#并行进程数
parser.add_argument("-p","--phe",help="phenotype", type=str)
parser.add_argument("-e","--exp",help="express", type=str,required=True)
parser.add_argument("-s","--type",help="Type(cds or prm)", type=str,required=True)
parser.add_argument("-o","--output",help="output file name", type=str,required=True)
args = parser.parse_args()
out_file = args.output
inputfile = args.inputfile
thread_num = args.threat
exp = args.exp
phe = args.phe
samples,genelist = get_sample(inputfile)
write_genefasta(samples,genelist)
folder1 = os.path.exists("var_file")
folder2 = os.path.exists("stat_file")
folder3 = os.path.exists("aln_file")
try:
    if not folder3:
        os.makedirs("aln_file")  
except:
    print("aln_file have exist")
try:
    if not folder1:
        os.makedirs("var_file") 
except:
    print("var_file  exist")
try :
    if not folder2:
        os.makedirs("stat_file") 
except:
    print("stat_file  exist")   
with open(out_file+".result","w") as f:
    f.write("#Geneid\tPOS\tTYPE\tREF\tALT\t")
    for id in sample_list:
        f.write(id+"\t")
    f.write("Valid_ref_count\tValid_alt_count\tTtest_pvalue\n") 
type = args.type    
if type == "cds" :
    person =  trait(exp,phe)
pool = Pool(processes=thread_num)
in_alns = pool.map(aln,genelist)
def filter(out_file):
    pvalues = []
    newlines = []
    with open(out_file+".result","r") as f1:
        for line in f1:
            lines = line.strip().split()
            geneid = lines[0]
            pos = lines[1]
            Type = lines[2]
            ref = lines[3]
            alt = lines[4]
            tValid_alt_count = lines[-1]
            Valid_ref_count = lines[-2]
            pvalue = lines[-3]
            newline = geneid+"\t"+pos+"\t"+Type+"\t"+ref+"\t"+alt+"\t"+tValid_alt_count+"\t"+Valid_ref_count+"\t"+pvalue
            newlines.append(newline)
            pvalues.append(pvalue)
    pvalues.sort()
    posnum = len(pvalues)
    min = pvalues[int(posnum*0.01)]  
    print(min)    
    with open(out_file+".out","w") as f2:
        f2.write("#Geneid\tPOS\tTYPE\tREF\tALT\tValid_ref_count\tValid_alt_count\tTtest_pvalue\n")
        for line in newlines:
            pvalue = line.split()[-1] 
            if pvalue <= min:
                f2.write(line+"\n")
filter(out_file)










"""
1 - pull out individuals from across multiple tsv files
2 - reduce each file to thinned data (first and last value <difffac)
3 - if asked, replace NAs with priors

DEPENDENCIES

Python 2.7

USAGE

python pull_thin_tsv.0.5.1.py <config_file: default = pt.cfg> 

If no config_file provided, program looks for pt.cfg in same folder

FILES TO INCLUDE IN SAME FOLDER


OUTPUT




1 - pull out individuals from tsv files
2 - reduce each file to thinned data (first and last value <difffac)
3 - if asked, replace NAs with priors
v.0.1 2/22/11
v.0.2 2/26/11 - bug fix
v.0.3 3/2/11 - bug fix - was incorrectly skipping NA2prior if indicated "all" to pull
v.0.4 3/11/11 - Feature upgrade - On X chrom, replace NA with prior, dependent on sex. If sex not provided, then assume all female.
v.0.4.1 3/18/11 - If user indicates "All" in any case for indivs in pt.cfg file, then read correctly.
v.0.5 - Improve speed performance of thin
v.0.5.1 - bug fix. Need to take abs value of diffs in thin
v.0.5.2 - bug fix. Will correctly look for sex_all parameter in pt.cfg
v.0.6 - allow files from only 1 msg run; check for thinned files; convert tsv to r-qtl csv format
v.0.6.1 - fixed assignment of female genotypes and probabilities
v.0.6.3 - fixed bug that incorrectly assigned states 1 and 2 to bc, changed to BB, AB
v.0.6.4 - fixed 3 bugs - 1 - used chrom as local variable in f2_reduce_prob_slots 
                        2 - failed to convert str to int for sex
                        3 - saved incorrect # of inds for f2 if sex was not indicated for all indivs in intrasex.tsv were. Now assumes missing sex are female.
v.0.6.5 - fixed bug that did not recognize absence of p1 file and converted both bc and f2 for f2
v.0.6.6 - 1 - fixed bugs that caused incorrect calling of female f2 X chromosome genotypes
          2 - sort all files by individual name (natural sort)
v.0.6.7 - bug fix - was replacing indiv names that contained "NA" with prior (e.g.indivA2_TNATTG)
v.0.6.8 - bug fix - was testing sex as logical 0/1, rather than numeric.
v.0.6.9 - changed scoring of males in f2 crosses to folow scoring in bc. Examine only whether par1 > or < par2
v.0.6.10 - bug fix - open files in Universal new line mode ('rU', which prevents error and failure to read phenofile). Unfortunately, code is buried in try statements, which hides this error from view. Should revise code to discard try statements.
v.0.6.11 - bug fix - not correctly recognizing sex in sexfile and phenofile
"""

import sys
import getopt
import ConfigParser, os
import csv
import glob
import time
from operator import delitem 


#default values
difffac = 0.01
chroms = "all"


def main(argv=None):
        if argv is None:
#		print "hello"
                config_file = "pt.cfg"
        else:
                config_file = argv
        #get parameters from config file
        config = ConfigParser.SafeConfigParser(allow_no_value=True)
        config.read(config_file)
        filePar2 = config.get('Common','filepar2')
        sorted_filePar2=[]
        for i in glob.iglob("*par2*sorted"):
                sorted_filePar2.append(i[:-7])
        if filePar2 in sorted_filePar2:
                print "%s has been pre-sorted" %(filePar2)
                pass
        if filePar2 not in sorted_filePar2:
                print "Sorting %s" %(filePar2)
                sort_file(filePar2,'\t')
        filePar2 = filePar2 + ".sorted"
        global cross
        try:
                cross = config.get('Common','cross')
                cross = cross.lower()
                print "Cross type = %s" %(cross)
        except:
                print "No cross type indicated. I will assume a backcross"
                cross = 'bc'
                #filePar1 = None

        if cross == 'f2':
        
                try:
                        filePar1 = config.get('Common','filepar1')
                        sorted_filePar1=[]
                        for i in glob.iglob("*par1*sorted"):
                                sorted_filePar1.append(i[:-7])
                        if filePar1 in sorted_filePar1:
                                print "%s has been pre-sorted" %(filePar1)
                                pass
                        if filePar1 not in sorted_filePar1:
                                print "Sorting %s" %(filePar1)
                                sort_file(filePar1,'\t')
                        filePar1 = filePar1 + ".sorted"
                except:
                        print "No Par1 file processed for f2 cross. Run aborted."
                        sys.exit()
                        

        num_inds = count_lines(filePar2)

        try:
                indivs = config.get('Common','indivs')

                if ',' in indivs:
                        indivs = indivs.replace(" ","").split(',')
                elif 'all' in indivs.lower():
                        indivs = "All"
                else:
                        indivs = indivs.split()
        except:
                indivs = "All"

        print "Individuals = %s" %(indivs)

        try:
                sex_all = config.get('Common','sex_all')
                if 'f' in sex_all.lower() or '0' in sex_all:
                        sex_all = '0'
                        sex = sex_all2sex(num_inds,sex_all)
                        print "All individuals are female"
                elif 'm' in sex_all.lower() or '1' in sex_all:
                        sex_all = '1'
                        sex = sex_all2sex(num_inds,sex_all)
                        print "All individuals are male"
        except:
                try:
                        phenofile = config.get('Common','phenofile')
                        #if supplied, get sex of each individual (0 = F, 1 = M)
                        
                        sorted_phenofile=[]
                        for i in glob.iglob("phenofile.sorted"):
                                sorted_phenofile.append(i[:-7])
                        if phenofile in sorted_phenofile:
                                print "%s has been pre-sorted" %(phenofile)
                                pass
                        if phenofile not in sorted_phenofile:             
                                sort_file(phenofile,'\t')

                        sex = get_sex_phenofile(phenofile + ".sorted")

                except:
                        try:
                                sexfile = config.get('Common','sexfile')
                                
                                sorted_sexfile=[]
                                for i in glob.iglob("sexfile.sorted"):
                                        sorted_sexfile.append(i[:-7])
                                if sexfile in sorted_sexfile:
                                        print "%s has been pre-sorted" %(sexfile)
                                        pass
                                if sexfile not in sorted_sexfile:             
                                        sort_file(sexfile,'\t')
                                
                                
                                sex = get_sex_sexfile(sexfile + ".sorted")
                        except:
                                sex = []
                                get_num_inds = open(filePar2,'rU')
                                for line in get_num_inds:
                                        sex.append('0')

                                print "No file with sex supplied. I will assume all females"	

        global difffac
        try:
                difffac = config.get('Common','difffac')
                difffac = float(difffac)
                print "Data will be thinned with difffac = %s" %(difffac)
        except:
                difffac = 0.01

        try:
                chroms = config.get('Common','chroms')
                if ',' in chroms:
                        chroms = chroms.replace(" ","").split(',')
                else:
                        chroms = chroms.split()
                        print "Chromosomes to pull = %s" %(chroms)
        except:
                print "No chromosomes specified. I will pull all"
                chroms = "all"




        global auto_prior
        global X_prior
        try:
                auto_prior = config.get('Common','autosome_prior')
                auto_prior = float(auto_prior)
                try:	
                        X_prior = config.get('Common','X_prior')
                        X_prior = float(X_prior)
                except ConfigParser.NoOptionError:
                        print "No X prior requested. X prior set to autosomal prior."
                        X_prior = auto_prior		
        except ConfigParser.NoOptionError:
                print "No autosomal or X prior set"
                auto_prior = False
                X_prior = False

        print "Autosomal prior to replace NAs = %s" %(auto_prior)
        print "X chromosome prior to replace NAs = %s" %(X_prior)	



        ##Pull requested individuals from Par2
        ##check if existing pulled filePar2 == all filePar2
        pre_pulled_filePar2=[]
        for i in glob.iglob("*par2*pulled"):
                pre_pulled_filePar2.append(i[:-7])
#		print set(pre_pulled_filePar2), set(filePar2)
        if filePar2 in pre_pulled_filePar2:
                print "%s has been pre-pulled" %(filePar2)
                pass
        if filePar2 not in pre_pulled_filePar2:#if file not yet pulled
                if "All" not in indivs:#if user specified particular individuals
                        pull_idds(filePar2,indivs)
                elif "All" in indivs:
                        file_to_pulled(filePar2)#if request all, just paste file to new file with .pulled suffix
        pre_converted_filePar2=[]
        for i in glob.iglob("*par2*.pulled.converted"):
                pre_converted_filePar2.append(i[:-17])
#		print set(pre_pulled_filePar2), set(filePar2)
        if filePar2 in pre_converted_filePar2:
                print "%s has been pre-converted" %(filePar2)
                pass
        if filePar2 not in pre_converted_filePar2:#if file not yet pulled
                NA2prior(filePar2,chroms,sex)


        ##If Par1 present, then pull requested individuals from Par1
        ##check if existing pulled filePar1 == all filePar1
        if cross == 'f2':
                pre_pulled_filePar1=[]
                for i in glob.iglob("*par1*pulled"):
                        pre_pulled_filePar1.append(i[:-7])
        #		print set(pre_pulled_filePar2), set(filePar2)
                if filePar1 in pre_pulled_filePar1:
                        print "%s has been pre-pulled" %(filePar1)
                        pass
                if filePar1 not in pre_pulled_filePar1:#if file not yet pulled
                        if "All" not in indivs:#if user specified particular individuals
                                pull_idds(filePar1,indivs)
                        elif "All" in indivs:
                                file_to_pulled(filePar1)#if request all, just paste file to new file with .pulled suffix
                pre_converted_filePar1=[]
                for i in glob.iglob("*par1*.pulled.converted"):
                        pre_converted_filePar1.append(i[:-17])
        #		print set(pre_pulled_filePar1), set(filePar1)
                if filePar1 in pre_converted_filePar1:
                        print "%s has been pre-converted" %(filePar1)
                        pass
                if filePar1 not in pre_converted_filePar1:#if file not yet pulled
                        NA2prior(filePar1,chroms,sex)


        ##Thin files based on numbers in Par2
        thinned_filePar2=[]#list of thinned filePar2
        for i in glob.iglob("*par2*pulled.converted.thinned"):
                thinned_filePar2.append(i[:-25])

        if filePar2 in thinned_filePar2:
                print "Parent 2 has been thinned already."
        elif filePar2 not in thinned_filePar2:
                if difffac:
                        if auto_prior or X_prior:
                                col_del = thin(filePar2)#thin file with difffac, return col_del for potential use to thing Par1
                        else:
                                col_del = thin_NA(filePar2)#thin file with difffac, return col_del for potential use to thing Par1			
                else:
                        print "No difffac provided. Data not thinned."
        ##Then use index of thinned markers to thin Par1
        thinned_filePar1=[]#list of thinned filePar2
        for i in glob.iglob("*par1*pulled.converted.thinned"):
                thinned_filePar1.append(i[:-25])
        renamed_filePar1=[]
        if cross == 'f2':
                if filePar1 in thinned_filePar1:
                        print "Parent 1 has been thinned already."
                elif filePar1 not in thinned_filePar1:#if not yet thinned
                        if difffac:
                                thin_Par1(filePar1,col_del)#thin file with difffac, return col_del for potential use to thing Par1

        ##Convert probs for f2
        if cross == 'f2':
                f2_filePar1 = []
                for i in glob.iglob("*par1*pulled.converted.thinned.f2_rqtl"):
                        f2_filePar1.append(i[:-33])

                if filePar1 in f2_filePar1:
                        print "Prob files have already been converted for f2 cross."
                else:
                        print "Converting X chrom data for f2 cross."
                        f2_reduce_prob_slots(filePar2,filePar1,sex)


        ##convert tsv to csv format, make hard calls
        if cross == 'bc':
                print "Creating csv file for bc"
                tsv2csv_bc(filePar2,sex)
        else:
                print "Creating csv file for f2 cross"
                tsv2csv_f2(filePar2,filePar1,sex)

# ---------------------------------------------------------
# natsort.py: Natural string sorting.
# ---------------------------------------------------------

# By Seo Sanghyeon.  Some changes by Connelly Barnes.

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp                
                
#################################


def sort_file(infile,delim):
        original_file = csv.reader(open(infile,'rU'),delimiter=delim)
        header = original_file.next()
        data_dict = dict([(line[0],line[1:]) for line in original_file])
        #sort by individual name

        sorted_dict_keys=natsorted(data_dict.keys())
        #write file
        outfile_name = infile + ".sorted"
        out_file = csv.writer(open(outfile_name,"w"),delimiter=delim)
        out_file.writerow(header)
        for x in sorted_dict_keys:
                row = data_dict[x]
                row.insert(0,x)
                out_file.writerow(row)
#        return outfile_name
                
                
def f2_reduce_prob_slots(par2,par1,sex):
        '''
        FOR FEMALES, CONVERT PAR2 PROB SLOTS TO PAR12 PROB
        THAT IS, PUT PROBPAR1 + PROBPAR2 IN NEWPROBPAR1 (homozygous)
        THEN, PUT PROBPAR12 (1-NEWPROBPAR1) IN PROBPAR2 (heterozygous)
        MALES STAY AS THEY ARE 
        '''
        thinned_file_par2 = csv.reader(open(par2 + ".pulled.converted.thinned", 'rU'),delimiter = '\t')
        thinned_file_par1 = csv.reader(open(par1 + ".pulled.converted.thinned", 'rU'),delimiter = '\t')

        converted_file_par2 = csv.writer(open(par2 + ".pulled.converted.thinned.f2_rqtl", "w"), delimiter = '\t')
        converted_file_par1 = csv.writer(open(par1 + ".pulled.converted.thinned.f2_rqtl", "w"), delimiter = '\t')

        header = thinned_file_par2.next()
        skip = thinned_file_par1.next()
        converted_file_par2.writerow(header)
        converted_file_par1.writerow(header)
        chromosomes = [i.split(":")[0] for i in header]
        ind = 0
        for row in thinned_file_par2:
                par1_line = thinned_file_par1.next()
                if sex[ind] == 1:#if male, calc prob_par1 = 1 - prob_par2
                        #converted_file_par2.writerow(row)
                        #converted_file_par1.writerow(par1_line)

                        #set up list to accept data
                        column = 0
                        converted_row_data_par2 = []
                        converted_row_data_par1 = []
                        for par2_datum in row:
                                par1_datum = par1_line[column]
                                #if chromosomes[column] != "X":
                                        #converted_row_data_par2.append(par2_datum)
                                        #converted_row_data_par1.append(par1_datum)
                                #elif chromosomes[column] == "X":
                                        #prob_par1 = 1 - float(par2_datum)
                                        #converted_row_data_par2.append(par2_datum)#het prob goes in par2 slot
                                        #converted_row_data_par1.append(par1_datum)#homo prob goes in par1 slot
                                converted_row_data_par2.append(par2_datum)
                                converted_row_data_par1.append(par1_datum)
                                column += 1

                        converted_file_par2.writerow(converted_row_data_par2)
                        converted_file_par1.writerow(converted_row_data_par1)

                       

                else:#if female, or missing (assumed female), more complicated
                        '''
			FOR FEMALES, CONVERT PAR2 PROB SLOTS TO PAR12 PROB
			THAT IS, PUT PROBPAR1 + PROBPAR2 IN NEWPROBPAR1 (homozygous)
                        THEN, PUT PROBPAR12 (1-NEWPROBPAR1) IN PROBPAR2 (heterozygous)
			'''
                        #set up list to accept data
                        column = 0
                        converted_row_data_par2 = []
                        converted_row_data_par1 = []
                        for par2_datum in row:
                                par1_datum = par1_line[column]
                                if column == 0:#if individual name
                                        converted_row_data_par2.append(par2_datum)
                                        converted_row_data_par1.append(par1_datum)
                                elif chromosomes[column] != "X":
                                        converted_row_data_par2.append(par2_datum)
                                        converted_row_data_par1.append(par1_datum)
                                elif chromosomes[column] == "X":
                                        prob_par1 = float(par2_datum) + float(par1_datum)
                                        prob_par12 = 1 - prob_par1
                                        converted_row_data_par2.append(prob_par12)#het prob goes in par2 slot
                                        converted_row_data_par1.append(prob_par1)#homo prob goes in par1 slot

                                column += 1

                        converted_file_par2.writerow(converted_row_data_par2)
                        converted_file_par1.writerow(converted_row_data_par1)

                ind += 1

def tsv2csv_bc(par2,sex):
        par2_thinned = par2 + ".pulled.converted.thinned"
        thinned_file = csv.reader(open(par2_thinned, 'rU'),delimiter = '\t')
        csv_out = csv.writer(open(par2 + ".csv", "w"), delimiter = ',')
        header = thinned_file.next()
        chroms = [i.split(":")[0] for i in header]
        header[0] = "id"
        csv_out.writerow(header)
        csv_out.writerow(chroms)
        csv_out.writerow("")

        x = 0#track # individuals
        for row in thinned_file:
                ind = row[0]
                markers = row[1:]
                ind_sex = sex[x]
                converted_markers = []
                converted_markers.append(ind)
                if ind_sex == '0':#if female, simple, genotypes = BA, BB
                        for marker in markers:
                                if marker == str(auto_prior):
                                        converted_markers.append("-")
                                elif float(marker) > 0.5:
                                        converted_markers.append("BB")#BB
                                else:
                                        converted_markers.append("BA")#AB
                else: #male
                        y = 1#track marker names
                        for marker in markers:
                                if chroms[y] != "X":
                                        if marker == str(auto_prior):
                                                converted_markers.append("-")
                                        elif float(marker) > 0.5:
                                                converted_markers.append("BB")#BB
                                        else:
                                                converted_markers.append("BA")#AB

                                elif chroms[y] == "X":
                                        if marker == str(X_prior):
                                                converted_markers.append("-")
                                        elif float(marker) > 0.5:
                                                converted_markers.append("BB")#BB
                                        else:
                                                converted_markers.append("AA")#AA			
                                y+=1
                x+=1

                csv_out.writerow(converted_markers)

def tsv2csv_f2(par2,par1,sex):
        par2_thinned = par2 + ".pulled.converted.thinned"#changed to sample from thinned, rather than f2_qtl data, to get female pgm calls correct
        par1_thinned = par1 + ".pulled.converted.thinned"
        thinned_file_par2 = csv.reader(open(par2_thinned, 'rU'),delimiter = '\t')
        thinned_file_par1 = csv.reader(open(par1_thinned, 'rU'),delimiter = '\t')
        csv_out = csv.writer(open(par2 + ".csv", "w"), delimiter = ',')
        header = thinned_file_par2.next()
        skip = thinned_file_par1.next()
        chroms_row = [i.split(":")[0] for i in header]
        chroms = chroms_row[1:]#this fixed the odd first X chrom marker being incorrectly scored v.0.6.6
        header[0] = "id"
        csv_out.writerow(header)
        csv_out.writerow(chroms_row)
        csv_out.writerow("")

        x = 0#track # individuals
        for row_p2 in thinned_file_par2:
                row_p1 = thinned_file_par1.next()
                ind = row_p2[0]
                markers_p2 = row_p2[1:]
                markers_p1 = row_p1[1:]
                ind_sex = sex[x]
                converted_markers = []
                converted_markers.append(ind)
                if ind_sex == 0:#if female, autosome genotypes = AA, AB, BB, X genotypes = AA (homozygote), AB

                        for z in xrange(len(markers_p2)):
                                genos = [float(markers_p1[z]), 1- float(markers_p2[z]) - float(markers_p1[z]), float(markers_p2[z])]#genotype probs for AA, AB, BB
                                geno = genos.index(max(genos))
                                if chroms[z] != "X":
                                        if genos == [auto_prior,auto_prior * 2,auto_prior]:
                                                converted_markers.append("-")
                                        elif geno == 0:
                                                converted_markers.append("AA")
                                        elif geno == 1:
                                                converted_markers.append("AB")
                                        elif geno == 2:
                                                converted_markers.append("BB")
                                elif chroms[z] == "X":
                                        if genos == [auto_prior,auto_prior * 2,auto_prior]:
                                                converted_markers.append("-")
                                        elif geno == 0:
                                                converted_markers.append("AA")#assign hard calls based on thinned data, but in Rqtl import pooled homozygote probs
                                        elif geno == 1:
                                                converted_markers.append("AB")
                                        elif geno == 2:
                                                converted_markers.append("BB")


                else: #male
                        for z in xrange(len(markers_p2)):
                                genos = [float(markers_p1[z]), 1- float(markers_p2[z]) - float(markers_p1[z]), float(markers_p2[z])]#genotype probs for AA, AB, BB
                                geno = genos.index(max(genos))
                                if chroms[z] != "X":
                                        if genos == [auto_prior,auto_prior * 2,auto_prior]:
                                                converted_markers.append("-")
                                        elif geno == 0:
                                                converted_markers.append("AA")#AA
                                        elif geno == 1:
                                                converted_markers.append("AB")#AB
                                        elif geno == 2:
                                                converted_markers.append("BB")#BB

                                elif chroms[z] == "X":
                                        #if genos == [auto_prior,auto_prior * 2,auto_prior]:
                                        #        converted_markers.append("-")
                                        if genos[0] == genos[2]:
                                                converted_markers.append("-")#if homozygotes have equal prob (for some reason, though this shouldn't happen)
                                        elif genos[0] > genos[2]:
                                                converted_markers.append("AA")#hemizygous A
                                        #elif geno == 1:
                                        #        converted_markers.append("-")#if male scored with het as highest prob, record as missing data
                                        elif genos[0] < genos[2]:
                                                converted_markers.append("BB")#hemizygous B

                x+=1

                csv_out.writerow(converted_markers)

def sex_all2sex(num_inds,sex_all):
        sex=[]
        for ind in range(num_inds):
                sex.append(sex_all)
        return sex

def get_sex_phenofile(phenofile):
        sex = []
        #get sex from phenofile
        open_file = csv.reader(open(phenofile,'rU'),delimiter='\t')
        header = open_file.next()
        sex_index = header.index('sex')
        for row in open_file:
                sex.append(int(row[sex_index]))
        replace_all(sex,'0',int('0'))
        replace_all(sex,'1',int('1'))
        replace_all(sex,'F',int('0'))
        replace_all(sex,'f',int('0'))
        replace_all(sex,'M',int('1'))
        replace_all(sex,'m',int('1'))
        return sex

def get_sex_sexfile(sexfile):
        sex = []
        #get sex from sexfile
        open_file = csv.reader(open(sexfile,'rU'),delimiter='\t')
        header = open_file.next()
        sex_index = header.index('sex')
        for row in open_file:
                sex.append(int(row[sex_index]))
        replace_all(sex,'0',int('0'))
        replace_all(sex,'1',int('1'))
        replace_all(sex,'F',int('0'))
        replace_all(sex,'f',int('0'))
        replace_all(sex,'M',int('1'))
        replace_all(sex,'m',int('1'))
        return sex

def with_index(seq):
        for i in xrange(len(seq)):
                yield i, seq[i]

def replace_all(seq, obj, replacement):
        for i, elem in with_index(seq):
                if elem == obj:
                        seq[i] = str(replacement)

def replace_all_in_array(array,obj,replacement):
        temp_array = []
        row_num = 0
        for row in array:
                row_num +=1
                if row_num==1:#ignore markers
                        temp_array.append(row)#just paste in temp_array
                else:
                        name = row[0]
                        data = row[1:]
                        replace_all(data,obj,replacement)
                        data.insert(0,name)
                        temp_array.append(data)
                        
        return temp_array

def pull_idds(sample_file,indivs):#filePar2,indivs):
        print "pulling requested individuals"
        open_file = open(sample_file,'rU')
        out_file = open(sample_file + ".pulled",'w')
        positions = open_file.readline()
        out_file.write(positions)
        for indiv in indivs:
                print indiv
                open_file = open(sample_file,'rU')
                for line in open_file:
                        if indiv in line:
                                ##CHECK FOR NOT ALL NA IN ROW
                                out_file.write(line)
        open_file.close()
        out_file.close()

def file_to_pulled(sample_file):
        open_file = csv.reader(open(sample_file,'rU'),delimiter='\t')
        print_file = csv.writer(open(sample_file + ".pulled",'w'),delimiter='\t')
        for line in open_file:
                print_file.writerow(line)

def NA2prior(sample_file,chroms,sex):
        #pull requested chromosomes
        #And convert NAs to priors
        print "Pulling requested chromosomes"
        temp_data=[]
        col_del=[]
        num_loci = 0
        num_ids = 0
        open_file = csv.reader(open(sample_file + ".pulled",'rU'),delimiter='\t')
        ##get data into lists
        for row in open_file:
                temp_data.append(row)
        if len(temp_data) > 1:#if any individuals in array
                ##find columns for deletions
                num_loci = len(temp_data[0])-1##note that loci start at pos 1
                num_ids = len(temp_data)


                ###TEST THIS SECTION
                #select only requested chromosomes
                chroms_lower=[x.lower() for x in chroms]
                if "all" not in chroms_lower:
                        selected_chrom_data = []
                        for row in temp_data:#save individuals from column 0
                                selected_chrom_data.append([row[0]])

                        for column in temp_data[0]:#takes column names from row 0
                                if column.split(":")[0] in chroms:#if the marker is from a requested chrom
                                        column_num = temp_data[0].index(column)#save index of marker
                                        #paste data from selected column into each row
                                        row_num = 0
                                        for row in temp_data:#for each row of data
                                                selected_chrom_data[row_num].append(row[column_num])#append data
                                                row_num +=1
                                        if column_num in range(1000,10000000,1000):
                                                print "Marker %s of %s pulled" %(column_num,num_loci)

                        temp_data = selected_chrom_data
                print "converting NAs to priors"

                #Replace NAs with priors, according to cross type, chromosome & sex
                if "bc" in cross:
                        if auto_prior:
                                
                                temp_data = replace_all_in_array(temp_data,"NA",auto_prior)
                        else:
                                print "No autosome prior provided for your backcross. NAs not replaced."
                elif "f2" in cross:
                        if X_prior and sex:
                                #replace male autosome with autosome_prior, X with 2*X_prior
                                #replace female autosome with autosome_prior, X with X+_prior
                                male_X_prior = 2 * X_prior
                                female_X_prior = X_prior
                                column_num=0#keep column index
                                for column in temp_data[0]:#takes column names from row 0
                                        if column_num >0:#ignore column of names
                                                if column.split(":")[0] in "X":#if marker on X
                                                        row_num = 0
                                                        for row in temp_data:#for each row of data
                                                                if row_num !=0:#ignore row 0, these are headers
                                                                        if "NA" in row[column_num]:
                                                                                if sex[row_num-1] == '1':#if male:Note, sex is array of sexes only, not header, thus reduce index by one relative to other arrays
                                                                                        temp_data[row_num][column_num] = male_X_prior
                                                                                else:#if female: Asume female if no sex parameters given
                                                                                        temp_data[row_num][column_num] = female_X_prior
        
                                                                row_num +=1
                                                column_num+=1

                                #Replace all autosomal NA's with autosome_prior
                                ##USE code from above for bc
                                if auto_prior:
                                        temp_data = replace_all_in_array(temp_data,"NA",auto_prior)
                                else:
                                        print "No autosome prior provided for your backcross. NAs not replaced."
                        else:
                                print "You failed to provide either/both X prior and sex for your F2 cross. NAs not replaced."




        converted_file = csv.writer(open(sample_file +".pulled.converted",'w'),delimiter='\t')
        for line in range(num_ids):
                converted_file.writerow(temp_data[line])



def thin(sample_file):

        print "collecting Par2 data for thinning"
        temp_data=[]
        col_del=[]
        num_loci = 0
        num_ids = 0
        open_file = csv.reader(open(sample_file + ".pulled.converted",'rU'),delimiter='\t')
        thinned_file = csv.writer(open(sample_file +".pulled.converted.thinned",'w'),delimiter='\t')

        ##get data into lists


        #Note that marker positions are now saved to separate array
        #post probs (and ind names) are save to temp_data. Makes easy to slice with temp_data[:][requested_column]
        idds = []
        markers = []
        markers = open_file.next()
        markers = markers[1:]#marker[0] = ''
        individual = 0
        for row in open_file:
                individual +=1
                if individual in range(10,10000,10):
                        print "Individual %s" %individual
                idds.append(row[0])
                pps = row[1:]
                pps = [float(a) for a in pps]
                temp_data.append(pps)
        row =[]
        #idds = [a[0] for a in temp_data]
        #temp_data = [a[1:] for a in temp_data]
        #for ind in temp_data:
        #	temp_data.append([float(a) for a in ind])
        num_loci = len(markers)
        num_ids = len(idds)
        start_time = time.clock()
        if num_ids > 1:#if any individuals in array
                ##find columns for deletions
                #num_loci = len(temp_data[1])##note that loci start at pos 1
                #num_ids = len(temp_data)##note that ids start at pos 1
                #identify redundant columns for deletion, defined by difffac
                print "identifying redundant columns"
                for marker in range(0,num_loci-2):
                        if marker in range(1000,10000000,1000):
                                print "Marker %s of %s" %(marker,num_loci)
                        if markers[marker].split(":")[0] == markers[marker+2].split(":")[0]: ##if same chrom arm

                                #get column of data for x, x+1, x+2
                                y = [i[marker] for i in temp_data]
                                y_1 = [i[marker+1] for i in temp_data]
                                y_2 = [i[marker+2] for i in temp_data]

                                if y == y_1 == y_2:#if 3 adjacent columns are identical
                                        col_del.append(marker+1)
                                else:
                                        Y1_diff=[]
                                        Y2_diff=[]

                                        Y1_diff = [abs(s - t) for s,t in zip(y_1,y)]
                                        Y2_diff = [abs(s - t) for s,t in zip(y_2,y_1)]	

                                        if max(Y1_diff) < difffac and max(Y2_diff) < difffac:
                                                col_del.append(marker+1)
                stop_time = time.clock()
                total_time = stop_time - start_time
                print "Time = %s" %total_time
                print "Deleting redundant columns from Par2"

                col_del.reverse()#need to delete from end of each line
                reverse_time = time.clock() - stop_time
                #print "Reversing = %s" %reverse_time
                #delete marked columns
                for del_marker in col_del:		
                        delitem(markers,del_marker)
                del_data=[]
                individual = 0
                for row in temp_data:
                        individual +=1
                        if individual in range(10,10000,10):
                                print "Individual %s" %individual 
                        for del_marker in col_del:
                                delitem(row,del_marker)
                        del_data.append(row)

                delete_time = time.clock() - reverse_time
                print "Deletion = %s" %delete_time
                print "Printing thinned file for Par2"
                #write thinned data to file
                markers.insert(0,'')
                thinned_file.writerow(markers)
                for line in range(num_ids):
                        del_data[line].insert(0,idds[line])
                        thinned_file.writerow(del_data[line])

        #return col_del for potential use to delete columns in Par1		
        return col_del

def thin_NA(sample_file):

        print "collecting Par2 data for thinning"
        temp_data=[]
        col_del=[]
        num_loci = 0
        num_ids = 0
        open_file = csv.reader(open(sample_file + ".pulled.converted",'rU'),delimiter='\t')
        thinned_file = csv.writer(open(sample_file +".pulled.converted.thinned",'w'),delimiter='\t')

        ##get data into lists


        #Note that marker positions are now saved to separate array
        #post probs (and ind names) are save to temp_data. Makes easy to slice with temp_data[:][requested_column]
        idds = []
        markers = []
        markers = open_file.next()
        markers = markers[1:]#marker[0] = ''
        individual = 0
        for row in open_file:
                individual +=1
                if individual in range(10,10000,10):
                        print "Individual %s" %individual
                idds.append(row[0])
                pps = row[1:]
                #pps = [float(a) for a in pps]
                temp_data.append(pps)
        row =[]
        #idds = [a[0] for a in temp_data]
        #temp_data = [a[1:] for a in temp_data]
        #for ind in temp_data:
        #	temp_data.append([float(a) for a in ind])
        num_loci = len(markers)
        num_ids = len(idds)
        start_time = time.clock()
        if num_ids > 1:#if any individuals in array
                ##find columns for deletions
                #num_loci = len(temp_data[1])##note that loci start at pos 1
                #num_ids = len(temp_data)##note that ids start at pos 1
                #identify redundant columns for deletion, defined by difffac
                print "identifying redundant columns"
                for marker in range(0,num_loci-2):
                        if marker in range(1000,10000000,1000):
                                print "Marker %s of %s" %(marker,num_loci)
                        if markers[marker].split(":")[0] == markers[marker+2].split(":")[0]: ##if same chrom arm

                                #get column of data for x, x+1, x+2
                                y = [i[marker] for i in temp_data]
                                y_1 = [i[marker+1] for i in temp_data]
                                y_2 = [i[marker+2] for i in temp_data]

                                if y == y_1 == y_2:#if 3 adjacent columns are identical
                                        col_del.append(marker+1)
                                #else:
                                #	Y1_diff=[]
                                #	Y2_diff=[]
        #
        #				Y1_diff = [s - t for s,t in zip(y_1,y)]
        #				Y2_diff = [s - t for s,t in zip(y_2,y_1)]	
        #						
        #				if max(Y1_diff) < difffac and max(Y2_diff) < difffac:
        #					col_del.append(marker+1)
                stop_time = time.clock()
                total_time = stop_time - start_time
                print "Time = %s" %total_time
                print "Deleting redundant columns from Par2"

                col_del.reverse()#need to delete from end of each line
                reverse_time = time.clock() - stop_time
                #print "Reversing = %s" %reverse_time
                #delete marked columns
                for del_marker in col_del:		
                        delitem(markers,del_marker)
                del_data=[]
                individual = 0
                for row in temp_data:
                        individual +=1
                        if individual in range(10,10000,10):
                                print "Individual %s" %individual 
                        for del_marker in col_del:
                                delitem(row,del_marker)
                        del_data.append(row)

                delete_time = time.clock() - reverse_time
                print "Deletion = %s" %delete_time
                print "Printing thinned file for Par2"
                #write thinned data to file
                markers.insert(0,'')
                thinned_file.writerow(markers)
                for line in range(num_ids):
                        del_data[line].insert(0,idds[line])
                        thinned_file.writerow(del_data[line])

        #return col_del for potential use to delete columns in Par1		
        return col_del

def thin_Par1(sample_file,col_del):

        temp_data=[]
        num_loci = 0
        num_ids = 0
        open_file = csv.reader(open(sample_file + ".pulled.converted",'rU'),delimiter='\t')
        thinned_file = csv.writer(open(sample_file +".pulled.converted.thinned",'w'),delimiter='\t')

        print "thinning data Par1"
        idds = []
        markers = []
        markers = open_file.next()
        markers = markers[1:]#marker[0] = ''
        individual = 0
        for row in open_file:
                individual +=1
                if individual in range(10,10000,10):
                        print "Individual %s" %individual
                idds.append(row[0])
                pps = row[1:]
                pps = [float(a) for a in pps]
                temp_data.append(pps)
        row =[]
        #idds = [a[0] for a in temp_data]
        #temp_data = [a[1:] for a in temp_data]
        #for ind in temp_data:
        #	temp_data.append([float(a) for a in ind])
        num_loci = len(markers)
        num_ids = len(idds)	


        #Delete columns based on redundant columns identified in Par2 file and passed to this function
        print "Deleting redundant columns from Par1"
        for del_marker in col_del:		
                delitem(markers,del_marker)
        del_data=[]
        individual = 0
        for row in temp_data:
                individual +=1
                if individual in range(10,10000,10):
                        print "Individual %s" %individual 
                for del_marker in col_del:
                        delitem(row,del_marker)
                del_data.append(row)

        print "Printing thinned file for Par1"
        markers.insert(0,'')
        thinned_file.writerow(markers)
        for line in range(num_ids):
                del_data[line].insert(0,idds[line])
                thinned_file.writerow(del_data[line])


def combine_filePar2(filePar2):
#Find all chromosome names, put in list
        chrom_positions=[]
        for file_name in filePar2:
                open_file = csv.reader(open(file_name + ".pulled.thinned",'rU'),delimiter='\t')
                row = open_file.next()
                chrom_positions.append(row)
        chrom_positions = [item for sublist in chrom_positions for item in sublist]#make single list sublists
        chrom_positions = list(set(chrom_positions))#reduce to unique set of chrom_positions

#Make list for each chromosome, then append all chromosome positions for all filePar2 in single list
#Check on chromosome names, then make set
        chroms = []
        for item in chrom_positions:
                if item:
                        chroms.append(item.split(":")[0])

                #sort(set(list)) to get non-redundant set of chrom positions
        chroms = sorted(set(chroms))
        chroms = [i + ":" for i in chroms]#replace colon, so that only real chromosomes (and not values including just numbers -like chr 4 - are called later)
        all_positions = ['']#header to paste in first row of csv file
        header=[]
        for chrom in chroms:
                chrom_bp=[]
                for item in chrom_positions:
                        if chrom in item:
                                chrom_bp.append(int(item.split(":")[1]))
                chrom_bp=sorted(set(chrom_bp))

                for position in chrom_bp:
                        all_positions.append(chrom + ":" + str(position))

        num_rows=0
        last_row = 0
        final_array=[]
        final_array.append(all_positions)

        for file_name in filePar2:
                ind_data=[]
                open_file = csv.reader(open(file_name + ".pulled.thinned",'rU'),delimiter='\t')
                file_header = open_file.next()
                for row in open_file:
                        ind_data.append(row)
                        final_array.append([row[0]])
                        num_rows+=1
                column_num=0
                for position in all_positions:
                        row_num=0+last_row
                        if position == '':
                                pass
                        elif position in file_header:
                                list_position = file_header.index(position)
                                for row in ind_data:
                                        row_num+=1
                                        final_array[row_num].append(row[list_position])
                        else:
                                for row in ind_data:
                                        row_num+=1
                                        final_array[row_num].append("NA")
                        column_num+=1
                last_row = row_num
        if "-par1.tsv" in filePar2[0]:
                suffix = "-par1.tsv"
        elif "-par2.tsv" in filePar2[0]:
                suffix = "-par2.tsv"
        elif "-par1par2.tsv" in filePar2[0]:
                suffix = "-par1par2.tsv"
        print_file = csv.writer(open("combined-ancestry-probs"+suffix,'w'),delimiter='\t')

        for line in final_array:
                print_file.writerow(line)

def count_lines(count_file):
        filename = open(count_file,'rU')
        x = 0
        for line in filename:
                x+=1
        return x
        filename.close()

if __name__ == "__main__":
        sys.exit(main())
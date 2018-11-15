"""
1 - pull out individuals from across multiple tsv files
2 - reduce each file to thinned data (first and last value <difffac)
3 - if asked, replace NAs with priors

DEPENDENCIES

Python 2.7

USAGE

python pull_thin_tsv.py <config_file: default = pt.cfg> 

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
v.0.6.12 - bug fix - now properly recognizes .cfg file passed on command line
v.0.6.13 - bug fix - Crashed if no prior assigned. Now force assignment to backcross prior if no prior assigned by user
v.0.6.14 - bug fix - Was not assigning X chromosme prior correctly for males
v.0.7 - implement core computations (convert and thin) in numpy - 1 Nov 2016
v.0.7.1 - fixed error when chrom has one marker, fixed assignment of sexes and Xchr when no sex and Xchr specific, fixed checking for par1 pre-runs
v.0.7.2 - fixed sex assignment when supplying sex in phenofile or sexfile
v.0.8 - added options to thin without considering missing data and to specify desired number of markers per chromosome, also thinning calculation is somewhat faster - 23 Oct 2018
"""

import sys
import getopt
import ConfigParser, os
import csv
import glob
import time
import re
import numpy as np
from operator import delitem 
import warnings
import random

warnings.filterwarnings("ignore",message = "All-NaN slice encountered")


#default values
#difffac = 0.01
chroms = "all"
xchroms = "X"

#def main(argv=None):
def main(argv=None):
        if argv is None:
                argv = sys.argv
#       print "hello"

                if len(argv) == 1:
                        config_file = "pt.cfg"
                else:
                        config_file = argv[1]
#        print argv
#        print config_file

        #get parameters from config file
         
         
        config = ConfigParser.SafeConfigParser(allow_no_value=True)
        config.read(config_file)
        filePar2 = config.get('Common','filepar2')
        sorted_filePar2=[]
        for i in glob.iglob("*par*sorted"):
                sorted_filePar2.append(i[:-7])
        if filePar2 in sorted_filePar2:
                print "%s has been pre-sorted" %(filePar2)
                pass
        if filePar2 not in sorted_filePar2:
                print "Sorting %s" %(filePar2)
                sort_file(filePar2,'\t')
        filePar2 = filePar2 + ".sorted"
        global cross
        if config.has_option('Common','cross'):
                cross = config.get('Common','cross')
                cross = cross.lower()
                print "Cross type = %s" %(cross)
        else:
                print "No cross type indicated. I will assume a backcross"
                cross = 'bc'
                #filePar1 = None

        if cross == 'f2':
        
                if config.has_option('Common','filepar1'):
                        filePar1 = config.get('Common','filepar1')
                        sorted_filePar1=[]
                        for i in glob.iglob("*par*sorted"):
                                sorted_filePar1.append(i[:-7])
                        if filePar1 in sorted_filePar1:
                                print "%s has been pre-sorted" %(filePar1)
                                pass
                        if filePar1 not in sorted_filePar1:
                                print "Sorting %s" %(filePar1)
                                sort_file(filePar1,'\t')
                        filePar1 = filePar1 + ".sorted"
                else:
                        print "No Par1 file processed for f2 cross. Run aborted."
                        sys.exit()
                        


        if config.has_option('Common','indivs'):
                indivs = config.get('Common','indivs')

                if ',' in indivs:
                        indivs = indivs.replace(" ","").split(',')
                elif 'all' in indivs.lower():
                        indivs = "All"
                else:
                        indivs = indivs.split()
        else:
                indivs = "All"

        print "Individuals = %s" %(indivs)

        num_inds = count_lines(filePar2) - 1

        if config.has_option('Common','sex_all'):
                sex_all = config.get('Common','sex_all')
                if 'f' in sex_all.lower() or '0' in sex_all:
                        sex_all = '0'
                        sex = sex_all2sex(num_inds,sex_all)
                        print "All individuals are female"
                elif 'm' in sex_all.lower() or '1' in sex_all:
                        sex_all = '1'
                        sex = sex_all2sex(num_inds,sex_all)
                        print "All individuals are male"
        else:
                print "No global sex indicated"

        if config.has_option('Common','phenofile'):
                phenofile = config.get('Common','phenofile')
                #if supplied, get sex of each individual (0 = F, 1 = M)
                
                sorted_phenofile=[]
                for i in glob.iglob(phenofile + ".sorted"):
                        sorted_phenofile.append(i[:-7])
                if phenofile in sorted_phenofile:
                        print "%s has been pre-sorted" %(phenofile)
                        pass
                else:
                        print "Sorting %s" %(phenofile)
                        sort_file(phenofile,'\t')
                ##Pull requested individuals from phenofile and/or sexfile
                pre_pulled_phenofile=[]
                for i in glob.iglob(phenofile + ".sorted.pulled"):
                        pre_pulled_phenofile.append(i[:-14])
                if phenofile in pre_pulled_phenofile:
                        print "%s has been pre-pulled" %(phenofile)
                        pass
                else:#if file not yet pulled
                        print "Pulling requested individuals from phenofile."
                        if "All" not in indivs:#if user specified particular individuals
                                pull_idds(phenofile + '.sorted',indivs)
                        elif "All" in indivs:
                                file_to_pulled(phenofile + '.sorted')#if request all, just paste file to new file with .pulled suffix       

                if config.has_option('Common','sex_all'):
                        pass 
                else:
                        sex = get_sex_phenofile(phenofile + ".sorted.pulled")

        elif config.has_option('Common','sexfile'):
                sexfile = config.get('Common','sexfile')
                sorted_sexfile=[]
                for i in glob.iglob(sexfile + ".sorted"):
                        sorted_sexfile.append(i[:-7])
                if sexfile in sorted_sexfile:
                        print "%s has been pre-sorted" %(sexfile)
                        pass
                else:
                        print "Sorting sexfile"
                        sort_file(sexfile,'\t')
                ##Pull requested individuals from phenofile and/or sexfile
                pre_pulled_sexfile=[]
                for i in glob.iglob(sexfile + ".sorted.pulled"):
                        pre_pulled_sexfile.append(i[:-14])
                if sexfile in pre_pulled_sexfile:
                        print "%s has been pre-pulled" %(sexfile)
                        pass
                else:#if file not yet pulled
                        print "Pulling requested individuals from sexfile."
                        if "All" not in indivs:#if user specified particular individuals
                                pull_idds(sexfile + '.sorted',indivs)
                        elif "All" in indivs:
                                file_to_pulled(sexfile + '.sorted')#if request all, just paste file to new file with .pulled suffix       
                if config.has_option('Common','sex_all'):
                        pass
                else:
                        sex = get_sex_phenofile(sexfile + ".sorted.pulled")
        else:
#                sex = []
                sex_all = '0'
                sex = sex_all2sex(num_inds,sex_all)
                print "No file with sex supplied. I will assume all females"                

        if config.has_option('Common','difffac') and config.has_option('Common','numMarkers'):
            print "Error, both difffac and numMarkers defined"
            sys.exit()
        
        global difffac
        if config.has_option('Common','difffac'):
                difffac = config.get('Common','difffac')
                difffac = float(difffac)
                numMarkers = None
                print "Data will be thinned with difffac = %s" %(difffac)
        elif config.has_option('Common','numMarkers'):
                difffac = None
                numMarkers = config.get('Common','numMarkers')
                numMarkers = int(numMarkers)
                print "Data will be thinned to %d markers per chromosome" %(numMarkers)
        else:
                difffac = 0.01
                numMarkers = None
                print "Data will be thinned with difffac = %s" %(difffac)
        
        if config.has_option('Common','ignoreNan'):
            tempHolder = config.get('Common','ignoreNan').lower()
            if tempHolder == "true" or tempHolder == "t":
                ignoreNans = True
                print "Ignoring missing data when thinning"
            elif tempHolder == "false" or tempHolder == "f":
                ignoreNans = False
            else: 
                ignoreNans = False
                print "Unrecognized 'ignoreNan' option. Defaulting to False"
        else:
            ignoreNans = False
                

        if config.has_option('Common','xchroms'):
                xchroms = config.get('Common','xchroms')
                if ',' in xchroms:
                        xchroms = xchroms.replace(" ","").split(',')
                else:
                        xchroms = xchroms.split()
                        print "Specified X chromosomes = %s" %(xchroms)
        else:
                print "No X chromosome specified. All chroms treated as autosomes"
                xchroms = "none"


        if config.has_option('Common','chroms'):
                chroms = config.get('Common','chroms')
                if ',' in chroms:
                        chroms = chroms.replace(" ","").split(',')
                else:
                        chroms = chroms.split()
                        print "Chromosomes to pull = %s" %(chroms)
        else:
                print "No chromosomes specified. I will pull all"
                chroms = "all"

        global auto_prior
        global X_prior
        if config.has_option('Common','autosome_prior'):
                auto_prior = config.get('Common','autosome_prior')
                auto_prior = float(auto_prior)
                if config.has_option('Common','X_prior'):   
                        X_prior = config.get('Common','X_prior')
                        X_prior = float(X_prior)
                else:
                        print "No X prior requested. X prior set to autosomal prior."
                        X_prior = auto_prior        
        else:
                print "No autosomal or X prior set. Both priors set to 0.5"
                auto_prior = 0.5
                X_prior = 0.5

        print "Autosomal prior to replace NAs = %s" %(auto_prior)
        print "X chromosome prior to replace NAs = %s" %(X_prior)   


        ##Pull requested individuals from Par2
        ##check if existing pulled filePar2 == all filePar2
        pre_pulled_filePar2=[]
        for i in glob.iglob("*par*pulled"):
                pre_pulled_filePar2.append(i[:-7])
#       print set(pre_pulled_filePar2), set(filePar2)
        if filePar2 in pre_pulled_filePar2:
                print "%s has been pre-pulled" %(filePar2)
                pass
        else:#if file not yet pulled
                if "All" not in indivs:#if user specified particular individuals
                        pull_idds(filePar2,indivs)
                elif "All" in indivs:
                        file_to_pulled(filePar2)#if request all, just paste file to new file with .pulled suffix
        filePar2 = filePar2 + ".pulled"
        
        
        
        pre_converted_filePar2=[]
        for i in glob.iglob("*par*.pulled.converted.thinned"):
                pre_converted_filePar2.append(i[:-18])
#       print set(pre_pulled_filePar2), set(filePar2)
        if filePar2 in pre_converted_filePar2:
                print "%s has been pre-converted and thinned" %(filePar2)
                pass
        else:#if file not yet converted and thinned
                print "Converting and thinning file"
                convert_and_thin(filePar2 ,chroms,sex,difffac,numMarkers,xchroms,ignoreNans)
                #NA2prior(filePar2,chroms,sex)
        filePar2 = filePar2 + ".converted.thinned"


        ##If Par1 present, then pull requested individuals from Par1
        ##check if existing pulled filePar1 == all filePar1
        if cross == 'f2':
                pre_pulled_filePar1=[]
                for i in glob.iglob("*par*pulled"):
                        pre_pulled_filePar1.append(i[:-7])
        #       print set(pre_pulled_filePar2), set(filePar2)
                if filePar1 in pre_pulled_filePar1:
                        print "%s has been pre-pulled" %(filePar1)
                        pass
                else:#if file not yet pulled
                        if "All" not in indivs:#if user specified particular individuals
                                pull_idds(filePar1,indivs)
                        elif "All" in indivs:
                                file_to_pulled(filePar1)#if request all, just paste file to new file with .pulled suffix
                filePar1 = filePar1 + ".pulled"

                pre_converted_filePar1=[]
                for i in glob.iglob("*par*.pulled.converted.thinned"):
                        pre_converted_filePar1.append(i[:-18])
        #       print set(pre_pulled_filePar1), set(filePar1)
                if filePar1 in pre_converted_filePar1:
                        print "%s has been pre-converted and thinned" %(filePar1)
                        pass
                else:#if file not yet pulled
                        print "Converting and thinning parent1"
                        convert_and_thin_Par1(filePar2,filePar1,sex,chroms,xchroms)
                        #convert_and_thin(filePar1,chroms,sex,difffac,xchroms)
                        #NA2prior(filePar1,chroms,sex)
                filePar1 = filePar1 + ".converted.thinned"


        ##Thin files based on numbers in Par2
        #thinned_filePar2=[]#list of thinned filePar2
        #for i in glob.iglob("*par*pulled.converted.thinned"):
                #thinned_filePar2.append(i[:-25])

        #if filePar2 in thinned_filePar2:
                #print "Parent 2 has been thinned already."
        #elif filePar2 not in thinned_filePar2:
                #if difffac:
                        #if auto_prior or X_prior:
                                #col_del = thin(filePar2)#thin file with difffac, return col_del for potential use to thin Par1
                        #else:
                                #col_del = thin_NA(filePar2)#thin file with difffac, return col_del for potential use to thin Par1          
                #else:
                        #print "No difffac provided. Data not thinned."
        ###Then use index of thinned markers to thin Par1
        #thinned_filePar1=[]#list of thinned filePar2
        #for i in glob.iglob("*par*pulled.converted.thinned"):
                #thinned_filePar1.append(i[:-25])
        #renamed_filePar1=[]
        #if cross == 'f2':
                #if filePar1 in thinned_filePar1:
                        #print "Parent 1 has been thinned already."
                #elif filePar1 not in thinned_filePar1:#if not yet thinned
                        #if difffac:
                                #thin_Par1(filePar1,col_del)#thin file with difffac, return col_del for potential use to thin Par1

        ##Convert probs for f2
        if cross == 'f2':
                f2_filePar1 = []
                for i in glob.iglob("*par*pulled.converted.thinned.f2_rqtl"):
                        f2_filePar1.append(i[:-33])

                if filePar1 in f2_filePar1:
                        print "Prob files have already been converted for f2 cross."
                else:
                        print "Converting X chrom data for f2 cross."
                        f2_reduce_prob_slots(filePar2,filePar1,sex)


        ##convert tsv to csv format, make hard calls
        sortedcsv = []
        for i in glob.iglob(filePar2 + ".csv"):
                sortedcsv.append(i[:-4])
        if filePar2 in sortedcsv:
                print "csv created alrady"
        else:
                if cross == 'bc':
                        print "Creating csv file for bc."
                        tsv2csv_bc(filePar2,sex)
                else:
                        print "Creating csv file for f2 cross."
                        tsv2csv_f2(filePar2,filePar1,sex)

        print "Done."
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
        
        
                
def f2_reduce_prob_slots(par2,par1,sex):
        '''
        FOR FEMALES, CONVERT PAR2 PROB SLOTS TO PAR12 PROB
        THAT IS, PUT PROBPAR1 + PROBPAR2 IN NEWPROBPAR1 (homozygous)
        THEN, PUT PROBPAR12 (1-NEWPROBPAR1) IN PROBPAR2 (heterozygous)
        MALES STAY AS THEY ARE
        Scale probabilities so sum to 1
        '''
        thinned_file_par2 = csv.reader(open(par2, 'rU'),delimiter = '\t')
        thinned_file_par1 = csv.reader(open(par1, 'rU'),delimiter = '\t')
        
        converted_file_par2 = csv.writer(open(par2 + ".f2_rqtl", "w"), delimiter = '\t')
        converted_file_par1 = csv.writer(open(par1 + ".f2_rqtl", "w"), delimiter = '\t')
        
        header = thinned_file_par2.next()
        skip = thinned_file_par1.next()
        converted_file_par2.writerow(header)
        converted_file_par1.writerow(header)
        chromosomes = [i.split(":")[0] for i in header]
        ind = 0
        for row in thinned_file_par2:
                par1_line = thinned_file_par1.next()
                if sex[ind] == '1':#if male, calc prob_par1 = 1 - prob_par2
                        #converted_file_par2.writerow(row)
                        #converted_file_par1.writerow(par1_line)
                        
                        #set up list to accept data
                        column = 0
                        converted_row_data_par2 = []
                        converted_row_data_par1 = []
                        for par2_datum in row:
                                par1_datum = par1_line[column]
                                if chromosomes[column] != "X":
                                        converted_row_data_par2.append(par2_datum)
                                        converted_row_data_par1.append(par1_datum)
                                elif chromosomes[column] == "X":
                                        prob_par1 = 1 - float(par2_datum)
                                        converted_row_data_par2.append(par2_datum)#het prob goes in par2 slot
                                        converted_row_data_par1.append(prob_par1)#homo prob goes in par1 slot
                                #converted_row_data_par2.append(par2_datum)
                                #converted_row_data_par1.append(par1_datum)
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
        par2_thinned = par2
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
        par2_thinned = par2#changed to sample from thinned, rather than f2_qtl data, to get female pgm calls correct
        par1_thinned = par1
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
                if ind_sex == '0':#if female, autosome genotypes = AA, AB, BB, X genotypes = AA (homozygote), AB
                        
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
        #thinned_file_par2.close()
        #thinned_file_par1.close()
        #csv_out.close()
        
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
                sex.append(str(row[sex_index]))
        replace_all(sex,0,'0')
        replace_all(sex,1,'1')
        replace_all(sex,'F','0')
        replace_all(sex,'f','0')
        replace_all(sex,'M','1')
        replace_all(sex,'m','1')
        return sex

def with_index(seq):
        for i in xrange(len(seq)):
                yield i, seq[i]

def replace_all(seq, obj, replacement):
        for i, elem in with_index(seq):
                if elem == obj:
                        seq[i] = str(replacement)

def pull_idds(sample_file,indivs):#filePar2,indivs):
        print "pulling requested individuals"
        open_file = open(sample_file,'rU')
        out_file = open(sample_file + ".pulled",'w')
        positions = open_file.readline()
        out_file.write(positions)
        for line in open_file:
                name = line.split()[0]
                for ind in indivs:
                        match = re.search(ind,name)
                        if match:
                                out_file.write(line)
                                break
        open_file.close()
        out_file.close()

def file_to_pulled(sample_file):
        open_file = csv.reader(open(sample_file,'rU'),delimiter='\t')
        print_file = csv.writer(open(sample_file + ".pulled",'w'),delimiter='\t')
        for line in open_file:
                print_file.writerow(line)
 
def convert_and_thin(sample_file,chroms,sex,difffac,numMarkers,xchroms,ignoreNans):
        #for testing
        #import numpy as np
        #sample_file="testpar1.tsv.sorted.pulled"
        #difffac = 0.1
        #chroms = ("1","X")
        #X_prior = 0.5
        #auto_prior = 0.25
        #sex = [1,1,1,0,0,0]
        #xchroms = 'X'
        
        #make numpy array of global sex
        #'sex' identifies sex of each individual
        sexarray = np.array(sex)#convert to numpy array for boolean search
        
        #get row and column #s
        with open(sample_file) as f:
                ncols = len(f.readline().split()) + 1
                nrows = sum(1 for _ in f)
        
        #read first row of *pulled tsv file to get chroms and positions
        print "reading markers"
        markers = np.genfromtxt(sample_file, max_rows=1, delimiter="\t",dtype='S')
        markers = np.delete(markers, 0, axis=0)
        
        #read first column
        print "reading individual names"
        inds = np.genfromtxt(sample_file,comments='#', delimiter="\t", 
                             usecols=0, usemask=False, loose=True, 
                             invalid_raise=True, max_rows=None, dtype='S')
        
        #read second column to end, skipping first row to get data
        
        #pp = np.genfromtxt(sample_file, skip_header=1,comments='#', delimiter="\t", missing_values="NA", 
                           #filling_values=np.nan, usecols=range(1,ncols), usemask=False, loose=True, 
                           #invalid_raise=True, max_rows=None)
        
        #get chromosomes from all markers
        ch = []
        m = np.char.rsplit(markers, sep=":", maxsplit=None)
        ch = np.array([x[0] for x in m]) #ch is list of chromosomes for each marker

        #IF CHROMS == ALL, LIST ALL CHROMS IN chroms
        if "all" in chroms:
                chroms = np.unique(ch)

        #PULL OUT ARRAY FOR EACH CHROMOSOME
        for chr in chroms:
                print "Thinning chromosome %s" %(chr)
                #chr_columns = np.char.find(markers, chr + ":") == 0
                chr_columns = ch == chr
                #SUBSET TOTAL ARRAY BASED ON INDICES IN CHR_COLUMNS
        
                marker_subset = markers[chr_columns]
#                ppsubset = np.array(pp[:,chr_columns])

                #GRAB CHROMOSOME FROM FILE
                chr_indices = np.where(ch==chr)
                chr_indices = np.squeeze(chr_indices)
                chr_indices = chr_indices + 1 #shift over one to allow for marker names
        
                pp = np.genfromtxt(sample_file, skip_header=1,comments='#', delimiter="\t", missing_values="NA", 
                                   filling_values=np.nan, usecols=chr_indices, usemask=False, loose=True, 
                                   invalid_raise=True, max_rows=None)

                try:
                        #markersthinned = np.empty((1,1),dtype = object)
                        #markersthinned[0,0] = marker_subset[0]

                        #if chrom contains only 1 marker, then place pp in ppthinned and skip thinning
                        np.shape(pp)[1]      #test if array has > 1 column, if this fails, will skip to except IndexError
        
                        #for each chrom, place first column
                        #ppthinned = np.empty((nrows, 1))
                        #ppthinned[:,0] = pp[:,0]
                
                        subset_ncols = np.shape(pp)[1]
                        #select appropriate intermediate columns
                        if not ignoreNans:
                            pp[np.isnan(pp)] = 2
                            difffac1 = np.max(np.abs(pp[:,0:-2] - pp[:,1:-1]),axis = 0)
                            difffac2 = np.max(np.abs(pp[:,0:-2] - pp[:,2:]),axis = 0)
                            pp[pp == 2] = np.nan
                            difffac1[difffac1 > 1] = 1 # Force the difffac to be 1 if comparing to an Nan
                            difffac2[difffac2 > 1] = 2 # Force the difffac to be 2 if comparing to an Nan (so after division, it's 1)
                        else:
                            difffac1 = np.nanmax(np.abs(pp[:,0:-2] - pp[:,1:-1]),axis = 0)
                            difffac2 = np.nanmax(np.abs(pp[:,0:-2] - pp[:,2:]),axis = 0)
                            
                            difffac1[np.isnan(difffac1)] = 0
                            difffac2[np.isnan(difffac2)] = 0
                        
                        if difffac != None:
                            difffac_list = np.maximum(difffac1,difffac2/2)
                            thinningLogical = difffac_list >= difffac
                        else: 
                            difffac_list = np.maximum(difffac1,difffac2/2)
                            sorted_list = np.sort(difffac_list)[::-1]
                            temp_difffac = sorted_list[numMarkers-3]
                            print "Thinning chr %s with difffac %f" % (chr,temp_difffac)
                            
                            if temp_difffac == 1:
                                print("Warning: the number of markers selected is less than the number of sites that switch from missing to not-missing data. You should set your number of markers higher. Thinning remaining markers randomly.")
                            
                            if np.sum(difffac_list == temp_difffac) > 1:
                                thinningLogical = difffac_list > temp_difffac
                                random.seed(42)
                                equalToDifffacIndices = numpy.where(difffac_list == temp_difffac)
                                selectedIndices = np.random.choice(equalToDifffacIndices,size=(numMarkers - np.sum(thinningLogical) - 2),replace=False)
                                thinningLogical[selectedIndices] = True
                            else:
                                thinningLogical = difffac_list >= temp_difffac
                        
                        thinningLogical = np.concatenate((True,thinningLogical,True),axis=None)
                        ppthinned = pp[:,thinningLogical]
                        markersthinned = np.array(marker_subset)[thinningLogical]
                        #ADD LAST COLUMN
                        #push = np.empty((nrows,1))#define array
                        
                        #if subset_ncols == 2: #only 2 columns, grab second
                        #        push[:,0] = pp[:,1]
                        #        markerpush = np.empty((1,1),dtype=object)
                        #        markerpush[0] = marker_subset[1]
                        #else:
                        #        push[:,0] = pp[:,col+2]#fill array with last column
                        #        markerpush = np.empty((1,1),dtype=object)
                        #        markerpush[0] = marker_subset[col+2]                                
                        #
                        #ppthinned = np.concatenate((ppthinned, push),axis=1)#concatenate column to ppthinned
                        #markersthinned = np.concatenate((markersthinned,markerpush),axis=0) 
                        
                except IndexError:
                        #array has one column, so skipped thinning and placed pp in ppthinned with correct shape. Proceed to conversion
                        ppthinned = np.repeat(pp, np.shape(markersthinned)[0],axis=0).reshape(-1,np.shape(markersthinned)[0])   
                        
                #convert NAs to priors
                if chr in xchroms: #X_prior
                        #print "sexarray=%s" %(sexarray)
                        malerows = sexarray=='1'
                        femalerows = sexarray=='0'
                        #print "malerows=%s" %(malerows)
                        #print "femalerows=%s" %(femalerows)
                        
                        malerows = np.repeat(malerows, np.shape(markersthinned)[0],axis=0).reshape(-1,np.shape(markersthinned)[0])
                        femalerows = np.repeat(femalerows, np.shape(markersthinned)[0],axis=0).reshape(-1,np.shape(markersthinned)[0])        
                        ppthinned[np.isnan(ppthinned) & malerows] = X_prior * 2
                        ppthinned[np.isnan(ppthinned) & femalerows] = X_prior
                else: #autosomal #auto_prior
                        ppthinned[np.isnan(ppthinned)] = auto_prior
        
        
                #print ppthinned
                #CONCATENATE ARRAYS FOR ALL CHROMOSOMES
                if chr == chroms[0]:
                        ppout = ppthinned
                        markers_out = markersthinned
                else:
                        ppout = np.concatenate((ppout,ppthinned),axis=1)
                        markers_out = np.concatenate((markers_out,markersthinned),axis=0)
        #    print ppout
        #    print markers_out
        
        #CONCATENATE MARKERS AND PP
        
        print "Joining marker names to genotype matrix"
        markedppout = np.concatenate((markers_out.reshape((1,markers_out.shape[0])),ppout),axis=0)
        inds = np.reshape(inds, (np.size(inds),1))
        outarray = np.concatenate((inds,markedppout),axis=1)
        #print outarray
        #WRITE FILE WITH PULLED CONVERTED THINNED DATA
        
        np.savetxt(sample_file +".converted.thinned", outarray, delimiter='\t', newline='\n',fmt = '%s')


def convert_and_thin_Par1(filePar2,filePar1,sex,chroms,xchroms):
        #grab first row of filePar2 and filePar2
        markersP2 = np.genfromtxt(filePar2, max_rows=1, delimiter="\t",dtype='S')
        markersP1 = np.genfromtxt(filePar1, max_rows=1, delimiter="\t",dtype='S')
        
        #compare and keep indices of filePar1 that match filePar2
        #grab_idx = np.searchsorted(markersP1,markersP2)
        
        grab_idx = []
        for i in range(np.size(markersP1)):
                if np.any(markersP1[i]==markersP2):
                        grab_idx.append(i)
        
        #read in sub_array based on grab_idx
        ppthinned = np.genfromtxt(filePar1, delimiter="\t", missing_values="NA", 
                           filling_values=np.nan, usecols=grab_idx, usemask=False, loose=True, 
                           invalid_raise=True, max_rows=None,dtype=object)
        #replace nan with prior
        #convert NAs to priors
        sexarray = np.array(sex)#convert to numpy array for boolean search
        
        ch = []
        m = np.char.rsplit(markersP2[1:], sep=":", maxsplit=None)
        for i in m:
                ch.append(i[0]) #ch is list of chromosomes for each marker
        ch = np.array(ch)
#        chroms = np.unique(ch) #ch is list of chromosomes
        
        #define male and female rows
        malerows = sexarray=='1'
        femalerows = sexarray=='0'
        #make arrays
        malerows = np.repeat(malerows, np.shape(markersP2)[0],axis=0).reshape(-1,np.shape(markersP2)[0])
        malerows = np.insert(malerows, 0, np.zeros(np.shape(markersP2)[0],dtype=bool)).reshape(-1,np.shape(markersP2)[0])#add column at top
        femalerows = np.repeat(femalerows, np.shape(markersP2)[0],axis=0).reshape(-1,np.shape(markersP2)[0])        
        femalerows = np.insert(femalerows, 0, np.zeros(np.shape(markersP2)[0],dtype=bool)).reshape(-1,np.shape(markersP2)[0])#add column at top
        
        #define X chrom columns
        ch = np.insert(ch,0,'') #insert blank value at beginning (genotype column)
        xch_idx = xchroms == ch
        #make array
        xchromcolumns = np.tile(xch_idx, np.shape(malerows)[0]).reshape(-1,np.shape(malerows)[1])
        
        ppthinned[(ppthinned == 'NA') & malerows & xchromcolumns] = X_prior * 2
        ppthinned[(ppthinned == 'NA') & femalerows & xchromcolumns] = X_prior
        
        #define all other NA as autosomal
        ppthinned[(ppthinned == 'NA')] = auto_prior
        
        #print out new thinned Par2
        np.savetxt(filePar1 +".converted.thinned", ppthinned, delimiter='\t', newline='\n',fmt = '%s')

def count_lines(count_file):
        filename = open(count_file,'rU')
        x = 0
        for line in filename:
                if line.strip():
                        x+=1
        filename.close()
        return x

if __name__ == "__main__":
        sys.exit(main())
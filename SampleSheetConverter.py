#! /usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys
import csv
import ConfigParser
import xml.etree.ElementTree as ET
import itertools
import shutil
import datetime

def create_dir(run, lane):
    flag_exit = 0

    if instrument.startswith('NS') or instrument.startswith('ST'):
        if not os.path.exists(os.path.join(run, 'Aligned_lane' + lane)):
            os.mkdir(os.path.join(run, 'Aligned_lane'+lane))
        else:
            print('ERROR: '+ run + '/Aligned_lane' + lane + ' already exists!', file=sys.stderr)
            flag_exit = 1

        if not os.path.exists(os.path.join(run, 'Unaligned_lane' + lane)):
            os.mkdir(os.path.join(run, 'Unaligned_lane' + lane))
            os.mkdir(os.path.join(run, 'Unaligned_lane' + lane, 'Undetermined_indices'))
        else:
            print('ERROR: ' + run + '/Unaligned_lane' + lane + ' already exists!', file=sys.stderr)
            flag_exit = 1
    else:
        if not os.path.exists(os.path.join(run, 'Aligned_lane' + lane)):
            os.mkdir(os.path.join(run, 'Aligned_lane' + lane))
        else:
            print('ERROR: ' +run + '/Aligned_lane' + lane + ' already exists!', file=sys.stderr)
            flag_exit = 1

        if not os.path.exists(os.path.join(run, 'Unaligned_lane' + lane)):
            os.mkdir(os.path.join(run, 'Unaligned_lane' + lane))
        else:
            print('ERROR: '+ run + '/Unaligned_lane' + lane + ' already exists!', file=sys.stderr)
            flag_exit = 1
    return flag_exit


def proc_sample_sheet(run, lane, iY, c_mode):
    flag_exit = 0
    sheet_orig = os.path.join(run, 'SampleSheetOriginal.csv')
    sheet_final = os.path.join(run, 'SampleSheet_lane' + lane + '.csv')
    bmode = set()
    flowcell_set = set()
    barcodes = dict()
    genomes = dict()
    with open(sheet_orig) as f_in:
        f_in_csv = csv.DictReader(f_in, dialect='excel')
        for f in f_in_csv:
            if f['Lane'] == lane:
                flowcell_set.add(f['FCID'])
                if len(f['Recipe']):
                    bmode.add(f['Recipe'].lower())
                else:
                    bmode.add(None)
                if f['SampleID'] not in barcodes.keys():
                    barcodes[f['SampleID']] = f['Index']
                    genomes[f['SampleID']] = f['SampleRef']
                else:
                    print('ERROR: '+ f['SampleID'] + ' is not unique!', file=sys.stderr)
                    flag_exit = 1

    if len(flowcell_set) > 1:
        print('ERROR: Different flowcells specified!', file=sys.stderr)
        flag_exit = 1
    if len(bmode) > 1:
        print('ERROR: Different barcoding modes (illumina, custom, no barcoding) used for the same lane!', file=sys.stderr)
        flag_exit = 1
    if not len(barcodes.values()) == len(set(barcodes.values())):
        print('ERROR: Barcodes are not unique!', file=sys.stderr)
        flag_exit = 1
    if len(set([len(bval) for bval in barcodes.values()])) > 1:
        print('ERROR: Barcodes are not of same length!', file=sys.stderr)
        flag_exit = 1
    fc = next(iter(flowcell_set))

    # Depending on the instrument type, generate SampleSheet.csv
    if instrument.startswith('NS') or instrument.startswith('ST') :
        with open(sheet_orig) as f_in, open(os.path.join(run, 'SampleSheet_lane' + lane + '.csv'), 'w') as samplesheet:
	    # Create SampleSheet special header
            f_out_cfg = ConfigParser.ConfigParser(allow_no_value = True)
            f_out_cfg.optionxform = str
            f_out_cfg.add_section('Header')
            f_out_cfg.set('Header', 'Investigator Name,GeneCore')
            f_out_cfg.set('Header', 'Project Name,' + fc)
            f_out_cfg.set('Header', 'Experiment Name,GeneCoreExp')
            f_out_cfg.set('Header', 'Date,'+ datetime.datetime.now().strftime("%m/%d/%Y"))
            f_out_cfg.set('Header', 'Workflow,GenerateFASTQ')
            f_out_cfg.add_section('Settings')
            f_out_cfg.add_section('Data')
            f_out_cfg.set('Data', ','.join(['SampleID', 'SampleName', 'index', 'index2']))

            f_in_csv = csv.DictReader(f_in, dialect='excel')
            for f in f_in_csv:
                if f['Lane'] == lane:
                    if f['Recipe'].lower() == "custom":
                        block = []
                        block.append([''.join(["lane", f['Lane'], ",", "lane", f['Lane'], ",", "", ",", ""])])
                        seen = []
                        for b in block:
                            if b in seen:
                                continue
                            seen.append(b)
                            f_out_cfg.set('Data', ''.join(b))
                    else:
                    	blen = max([len(bval) for bval in barcodes.values()])
                        if len(iY) == 2 and blen > int(c_mode[iY[0]]): #2nd condition to check that indeed dual index mode
                            f_out_cfg.set('Data', ','.join([f['SampleID'], f['SampleID'], f['Index'][:len(f['Index'])/2], f['Index'][len(f['Index'])/2:]]))
                        else:
                            f_out_cfg.set('Data', ','.join([f['SampleID'], f['SampleID'], f['Index'], ""]))

            # Writing Samplesheet.csv file
            f_out_cfg.write(samplesheet)

	# Copy SampleSheet.csv
        if instrument.startswith('NS') and lane == '1':
            shutil.copyfile(os.path.join(run, 'SampleSheet_lane' + lane + '.csv'), os.path.join(run, 'SampleSheet.csv'))
    else:
        with open(sheet_orig) as f_in, open(sheet_final, 'w') as f_out:
            fieldnames = ['FCID', 'Lane', 'SampleID', 'SampleRef', 'Index', 'Description', 'Control', 'Recipe', 'Operator', 'SampleProject']
            f_in_csv = csv.DictReader(f_in, dialect='excel')
            f_out_csv = csv.DictWriter(f_out, fieldnames=fieldnames)
            f_out_csv.writeheader()

            # Create SampleSheet with custom barcoding lanes removed
            for f in f_in_csv:
                if f['Lane'] == lane:
                    if f['Recipe'].lower() != "custom":
                        f_out_csv.writerow({'FCID': f['FCID'], 'Lane':f['Lane'], 'SampleID':f['SampleID'], 'SampleRef':f['SampleRef'], 'Index':f['Index'], 'Description':f['Description'], 'Control':f['Control'], 'Recipe':f['Recipe'].lower(), 'Operator':f['Operator'], 'SampleProject':f['SampleProject']})
                    else:
                        # We allow only one barcoding scheme per lane so the entire lane is custom barcoding
                        f_out_csv.writerow({'FCID':f['FCID'], 'Lane':f['Lane'], 'SampleID': "lane"+f['Lane'], 'SampleRef':f['SampleRef'], 'Index':'', 'Description':f['Description'], 'Control':f['Control'], 'Recipe':'', 'Operator':f['Operator'], 'SampleProject':f['SampleProject']})
                        break

    # Create directories for each samples
    if instrument.startswith('NS') or instrument.startswith('ST') :
        os.mkdir(os.path.join(run, 'Unaligned_lane' + lane, 'Project_' + fc))
    return fc, next(iter(bmode)), barcodes, genomes, flag_exit

def process_runinfo(run):
    tree = ET.parse(os.path.join(run, 'RunInfo.xml'))
    root = tree.getroot()
    # Instrument
    instr = root[0][1].text
    idx_mode = []
    cyc_mode = []
    for r in root.findall("./Run/Reads/"):
        idx_mode.append(r.get('IsIndexedRead'))
        cyc_mode.append(r.get('NumCycles'))
    return instr, idx_mode, cyc_mode

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def proc_base_mask(indN, indY, cyc, bmode, barcodes):
    flag_exit = 0
    blen = 0
    if (bmode is not None) and (bmode.upper() != 'CUSTOM'):
        blen = max([len(bval) for bval in barcodes.values()])
    if len(indN) > 1:
        pr = 'PE'
        if (len(indY) == 0) and (blen == 0):
            usebases = 'Y' + cyc[indN[0]] + ',' + 'Y' + cyc[indN[1]]
        elif len(indY) == 1: # Paired-end - 1 barcode
            usebases = 'Y' + cyc[indN[0]] + ','+ 'I'*blen + 'n'*(int(cyc[indY[0]])-blen) + ',' + 'Y' + cyc[indN[1]]
        elif len(indY) == 2: #Paired-end - dual barcodes
        	if blen <= int(cyc[indY[0]]): #Check if bmode is indeed dualbarcoded
        		usebases = 'Y' + cyc[indN[0]] + ',' + 'I'*(blen) + 'n'*(int(cyc[indY[0]])-blen) + ',' + 'n'*(int(cyc[indY[1]])) + ',' + 'Y' +cyc[indN[1]]
        	else:
        		usebases = 'Y' + cyc[indN[0]] + ',' + 'I'*(blen/2) + 'n'*(int(cyc[indY[0]])-blen/2) + ',' + 'I'*(blen/2) + 'n'*(int(cyc[indY[1]])-blen/2) + ',' + 'Y' +cyc[indN[1]]
        else:
		print('ERROR: Basemask error!', file=sys.stderr)
        	flag_exit = 1
        	
    else:
        pr = 'SE'
        if len(indY) == 0:
			usebases = 'Y' + cyc[indN[0]]
        elif (len(indY) == 1):
			if (blen == 0): # Single-end custom barcoded 
				usebases = 'Y' + cyc[indN[0]] + ',' + 'n' + cyc[indY[0]]
			else:
				usebases = 'Y' + cyc[indN[0]] + ',' + 'I'*blen + 'n'*(int(cyc[indY[0]])-blen)
        elif len(indY) == 2: #Single-end + dual barcodes
		if blen <= int(cyc[indY[0]]): #Check if bmode is indeed dualbarcoded
			usebases = 'Y' + cyc[indN[0]] + ',' + 'I'*(blen) + 'n'*(int(cyc[indY[0]])-blen) + ',' + 'n'*(int(cyc[indY[1]]))
		else:
            		usebases = 'Y' + cyc[indN[0]] + ',' + 'I'*(blen/2) + 'n'*(int(cyc[indY[0]])-blen/2) + ',' + 'I'*(blen/2) + 'n'*(int(cyc[indY[1]])-blen/2)
        else:
            print('ERROR: Basemask error!', file=sys.stderr)
            flag_exit = 1
    return pr, usebases, flag_exit

# Parse command line
parser = argparse.ArgumentParser(description='Genecore Illumina pipeline.')
parser.add_argument('-r', '--run', metavar='DIR', required=True, dest='runfolder', help='run folder')
parser.add_argument('-l', '--lane', metavar='1', required=True, dest='lane', help='lane number')
parser.add_argument('-c', '--config', metavar='SampleSheetConverter.cfg', required=True, dest='config', help='config file name')
args = parser.parse_args()

# Force RunVol mount
if os.system('stat {} >/dev/null'.format(args.runfolder)) != 0:
    sys.exit(1)

# Parse config file
cfg = ConfigParser.RawConfigParser()
cfg.read(args.config)

# Parsing RunInfo
print('Parsing RunInfo', file=sys.stdout)
instrument, index_mode, cycle_mode = process_runinfo(args.runfolder)
print('Instrument: ' + instrument + ', Index-mode: ' + '-'.join(str(x) for x in index_mode) +', Cycle-mode: ' + '-'.join(str(c) for c in cycle_mode), file=sys.stdout)
indexes_N = [i for i, x in enumerate(index_mode) if x == 'N']
indexes_Y = [i for i, x in enumerate(index_mode) if x == 'Y']

# Create Directories
print('Creating Aligned and Unaligned directories.', file=sys.stdout)
exit_dir = create_dir(args.runfolder, args.lane)
shutil.copyfile(args.config, os.path.join(args.runfolder, 'SampleSheet_lane' + args.lane + '.cfg'))

# Parsing sample sheet
print('Parsing sample sheet.', file=sys.stdout)
flowcell, barcode_mode, barcodes_ss, genomes_ss, exit_sheet = proc_sample_sheet(args.runfolder, args.lane, indexes_Y, cycle_mode)
print('Flowcell: ' + flowcell + ', Lane: ' + args.lane + ', Barcode-mode: ' + (barcode_mode if barcode_mode else 'None'), file=sys.stdout)

# Generating base_mask
print('Get base_mask info.', file=sys.stdout)
paired, use_bases, exit_basemask = proc_base_mask(indexes_N, indexes_Y, cycle_mode, barcode_mode, barcodes_ss)
print('Paired/Single-end: ' + paired + ', Basemask: ' + use_bases, file=sys.stdout)

# Check if genome is present
print('Genome checks.', file=sys.stdout)
exit_genome = 0
for g in set(genomes_ss.values()):
    if (g.lower() != 'unknown') and (not os.path.exists(cfg.get('dirs', 'genomes')+'/'+g)):
        print('ERROR: Genome ' + g +' is not present in '+ cfg.get('dirs', 'genomes'), file=sys.stdout)
        exit_genome = 1

# Check hamming distance
hamming = 'NA'
if (len(barcodes_ss) > 1):     
     if barcode_mode is not None:
	hamming = str(min([hamming_distance(e[0], e[1]) for e in itertools.combinations(barcodes_ss.values(), 2)]))

# Write infos in cfg file
cfg_new = ConfigParser.RawConfigParser()
cfg_new.add_section('PipelineInfo')
cfg_new.set('PipelineInfo', 'lane', args.lane)
cfg_new.set('PipelineInfo', 'flowcell', flowcell)
cfg_new.set('PipelineInfo', 'barcode_mode', (barcode_mode if barcode_mode else 'None'))
cfg_new.set('PipelineInfo', 'instrument', instrument)
cfg_new.set('PipelineInfo', 'paired', paired)
cfg_new.set('PipelineInfo', 'use_bases', use_bases)
cfg_new.set('PipelineInfo', 'genome', ",".join(set(genomes_ss.values())))
cfg_new.set('PipelineInfo', 'minimum_Hamming_distance', hamming)
with open(os.path.join(args.runfolder, 'SampleSheet_lane'+args.lane+'.cfg'), 'a') as configfile:
    cfg_new.write(configfile)

# Any errors?
if (exit_dir) or (exit_sheet) or (exit_basemask) or (exit_genome):
    #bc=open(os.path.join(run, 'error_cookie' + lane + '.txt'),"w")
    #bc.close()
    sys.exit(1)
else:
    print('All samples written.', file=sys.stdout)
    print('Completed without errors.', file=sys.stdout)

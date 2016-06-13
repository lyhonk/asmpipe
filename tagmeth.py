#!/usr/bin/env python

from __future__ import print_function

import re
import os
import os.path
import optparse
import sys
import pysam

# global regular expressions

readKeyRe = re.compile('(.*)\s[1-2](.*)')
readNumRe = re.compile('.*\s([1-2]).*')
cgPosRe = re.compile('(CG)')
tagCTRe = re.compile('F1/CT')
tagGARe = re.compile('F1/GA')

# sort bam file with pysam.sort function
# the bam file should be sorted by qname

def SortSam(inBam, outBam):
	pysam.sort("-n", inBam, outBam)

def TrimReadSeq(seq, cigar):
	trimedSeq = ''
	pos = 0
	for m in cigar :
		if(m[0] == 0) : # M
			trimedSeq += seq[pos : pos + m[1]]
			pos += m[1]
		elif(m[0] == 1): # I
			pos += m[1]
		elif(m[0] == 2): # D
			trimedSeq += ('N' * m[1])
		elif(m[0] == 3): # N
			trimedSeq += ('N' * m[1])
		elif(m[0] == 4): # S
			pos += m[1]
	return trimedSeq

def FormatPosList(posList):
	if(len(posList) == 0):
		return 'NA'
	return ','.join(map(str, posList))

def WriteTagMeth(tag, tagType, cvt, chrname, pos, tagLen, meth, unmeth, undt, outFile):
	outFile.write('%s\t%s\t%s\t%s\t%15ld\t%6d\t%s\t%s\t%s\t%d\t%d\t%d\n' % (tag, tagType, cvt, 
		chrname, pos, tagLen,
		FormatPosList(meth),
		FormatPosList(unmeth), 
		FormatPosList(undt),
		len(meth), len(unmeth), len(undt)))

# get tagmeth of single read

def TagMethOfSingleReads(dictSingle, dictRefSeq, bamFile, outFile):

	for read in dictSingle.itervalues():
		readSeq = TrimReadSeq(read.seq, read.cigar)
		tagLen = len(readSeq)
		readPos = read.pos
		tag = read.qname
		chrname = bamFile.getrname(read.rname)

		# get the corresponding part of reference sequence

		if(not (chrname in dictRefSeq)):
			print('\nwarning: chrom "' + chrname + '" was not found in the reference sequences')
			pass
		refSeq = dictRefSeq[chrname][readPos : readPos + tagLen]

		# CG positions on reference sequence

		cgPos = [ m.start() for m in cgPosRe.finditer(refSeq) ]

		# check CG methylation

		cvtCT = False
		cvtGA = False
		if(not read.has_tag('XB')):
			return
		tagXB = read.get_tag('XB')
		if(tagCTRe.match(tagXB)):
			cvtCT = True
			cvt = 'CT'
		elif(tagGARe.match(tagXB)):
			cvtGA = True
			cvt = 'GA'
		else:
			return

		meth = []
		unmeth = []
		undt = []
		
		if(cvtCT):
			for pos in cgPos :
				readBase = readSeq[pos]
				if (readBase == 'C'):
					meth += [ pos ]
				elif (readBase == 'T'):
					unmeth += [ pos ]
				else:
					undt += [ pos ]
		elif(cvtGA):
			for pos in cgPos :
				readBase = readSeq[pos + 1]
				if (readBase == 'G'):
					meth += [ pos ]
				elif (readBase == 'A'):
					unmeth += [ pos ]
				else:
					undt += [ pos ]

		if(not (len(meth) <= 0 and len(unmeth) <= 0 and len(undt) <= 0)):
			WriteTagMeth(tag.split()[0], 'S', cvt, chrname, readPos, tagLen, meth, unmeth, undt, outFile)
	

# get tagmeth of paired read

def TagMethOfPairedReads(dictPaired, dictRefSeq, bamFile, outFile):
	for pair in dictPaired.itervalues():
		if(not (pair[0] and pair[1])):
			return	

		chrname = bamFile.getrname(pair[0].rname)
		if(not (chrname in dictRefSeq)):
			print('\nwarning: chrom "' + chrname + '" was not found in the reference sequences')
			pass
		
		readSeq1 = TrimReadSeq(pair[0].seq, pair[0].cigar)
		readSeq2 = TrimReadSeq(pair[1].seq, pair[1].cigar)
		readSeqLen1 = len(readSeq1)
		readSeqLen2 = len(readSeq2)
 
		tagPos = pair[0].pos
		tagLen = pair[1].pos + readSeqLen2 - tagPos
		refSeq = dictRefSeq[chrname][tagPos : tagPos + tagLen]
		if(tagLen >= readSeqLen1 + readSeqLen2):
			tagSeq = readSeq1 + 'N' * (pair[1].pos - pair[0].pos - readSeqLen1) + readSeq2
		elif(tagLen > readSeqLen1 and tagLen < (readSeqLen1 + readSeqLen2)):
			tagSeq = readSeq1 + readSeq2[readSeqLen1 + readSeqLen2 - tagLen : ]
		else:
			tagSeq = readSeq1[:tagLen]

		cgPos = [ m.start() for m in cgPosRe.finditer(refSeq) ]		
		
		if(not (pair[0].has_tag('XB') and pair[1].has_tag('XB'))):
			return
		if(not (pair[0].get_tag('XB') == pair[1].get_tag('XB'))):
			return
		tagXB = pair[0].get_tag('XB')

		cvtCT = False
		cvtGA = False
		if(tagCTRe.match(tagXB)):
			cvtCT = True
			cvt = 'CT'
		elif(tagGARe.match(tagXB)):
			cvtGA = True
			cvt = 'GA'
		else:
			return

		meth = []
		unmeth = []
		undt = []
		if(cvtCT):
			for pos in cgPos :
				readBase = tagSeq[pos]
				if (readBase == 'C'):
					meth += [ pos ]
				elif (readBase == 'T'):
					unmeth += [ pos ]
				else:
					undt += [ pos ]
		elif(cvtGA):
			for pos in cgPos :
				readBase = tagSeq[pos + 1]
				if (readBase == 'G'):
					meth += [ pos ]
				elif (readBase == 'A'):
					unmeth += [ pos ]
				else:
					undt += [ pos ]

		if(not (len(meth) <= 0 and len(unmeth) <= 0 and len(undt) <= 0)):
			WriteTagMeth(pair[0].qname.split()[0], 'P', cvt, chrname, tagPos, tagLen, meth, unmeth, undt, outFile)

def main():

	# parse the command line options

	usage = 'usage: %prog [options] input.bam refseq.fa -o output.csv'
	parser = optparse.OptionParser(usage=usage, version='%prog 0.1.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input BAM file before data processing')
	parser.add_option('-q', '--quiet-mode', 
						action="store_true", dest="quiet", default=False,
						help='supress progress update information')

	(options, args) = parser.parse_args()
	if(len(args) != 2):
		parser.print_help()
		sys.exit(0)
	
	inputBamFileName = args[0]
	refSeqFileName = args[1]
	bamFileName = inputBamFileName
	baseFileName = os.path.splitext(os.path.basename(inputBamFileName))[0]
	outputFileName =  baseFileName + '.tagmeth.csv'
	logFileName = baseFileName + '.tagmeth.log'

	# sort the input bam file if -s option is set

	rmTemp = False
	if(options.sort):
		print('[*] Sorting by QNAME...')
		bamFileName = 'sorted.' + baseFileName
		SortSam(inputBamFileName, bamFileName)
		bamFileName += '.bam'
		rmTemp = True

	# load input files

	print('[*] Initializing...')

	if(not os.path.exists(bamFileName)):
		print('error: Failed to open file "', bamFileName, '"')
		sys.exit(-1)
	bamFile = pysam.AlignmentFile(bamFileName, "rb")
	
	if(not os.path.exists(refSeqFileName)):
		print('error: Reference sequence file "', refSeqFileName, '"', ' doest not exist.')
		sys.exit(-1)

	dictRefSeq = {}
	with open(refSeqFileName, 'r') as refSeqFile :
		chrname = ''
		seq = ''
		for line in refSeqFile:
			if(line[0] == '>'):
				
				# save current seq for current chr

				if(chrname != ''):
					dictRefSeq[chrname] = seq
				
				# new chrname & seq

				chrname = line[1:].strip()
				seq = ''
				print('    loading reference sequence: ' + chrname)
			else:
				seq += line.strip().upper()
		
		# write the last chr

		if(chrname != ''):
			dictRefSeq[chrname] = seq
	refSeqFile.close()

	# prepare output files

	if(options.outputfile):
		outputFileName = options.outputfile
	try:
		outFile = open(outputFileName, 'w')
	except IOError:
		print('error: write to output file failed!')
		sys.exit(-1)
	outFile.write('tagname\ttype\tcvt\tchr\tpos\tlen\tmethylated\tunmethylated\tundetermined\tmethycount\tunmethycount\tundtcount\n')

	# analyse algnments

	print('[*] Analyzing...')
	
	dictSingle = {}
	dictPaired = {}
	currentGroupKey = ''
	readCount = 0
	for read in bamFile.fetch(until_eof = True):
		try:
			groupKey = ''.join(readKeyRe.findall(read.qname)[0])
		except:
			groupKey = read.qname

		if(groupKey != currentGroupKey):
			currentGroupKey = groupKey

			# handle and write results

			TagMethOfPairedReads(dictPaired, dictRefSeq, bamFile, outFile)
			TagMethOfSingleReads(dictSingle, dictRefSeq, bamFile, outFile)
			dictPaired.clear()
			dictSingle.clear()

		# check if it is properly paired

		if(read.is_proper_pair):

			# paired reads

			chrpos = bamFile.getrname(read.rname).strip() + str(read.pos)
			chrposNext = bamFile.getrname(read.rnext).strip() + str(read.pnext)
			if(chrpos < chrposNext):
				readKey = groupKey + ':' + chrpos + ':' + chrposNext
			else:
				readKey = groupKey + ':' + chrposNext + ':' + chrpos

			if(read.is_reverse):
				readType = 1
			else:
				readType = 0
			
			if(readType == 0 or readType == 1):
				if(readKey in dictPaired):
					dictPaired[readKey][readType] = read
				else:
					dictPaired[readKey] = [None, None]
					dictPaired[readKey][readType] = read
		else:

			# single read

			try:
				readKey = groupKey + ':' + readNumRe.findall(read.qname)
			except:
				readKey = groupKey
			dictSingle[readKey] = read

		# progress

		readCount += 1
		if not options.quiet :
			sys.stdout.write('\r    read: #%ld' % (readCount))
			sys.stdout.flush()

	TagMethOfPairedReads(dictPaired, dictRefSeq, bamFile, outFile)
	TagMethOfSingleReads(dictSingle, dictRefSeq, bamFile, outFile)
	
	if not options.quiet :
		sys.stdout.write('\n')
	
	# release resources

	bamFile.close()
	outFile.close()
	
	if(rmTemp):
		os.unlink(bamFileName)

	print('[*] Complete')

if __name__ == '__main__':
	main()

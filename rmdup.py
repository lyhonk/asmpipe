#!/usr/bin/env python

from __future__ import print_function

import re
import sys
import os
import os.path
from optparse import OptionParser
import pysam

readKeyRe = re.compile('(.*)\s[1-2](.*)')
readNumRe = re.compile('.*\s([1-2]).*')

# upper and lower limit of the paired reads length

upperLimit = 2000
lowerLimit = 0

# [*] Log file Specification
# 
# Symbol	Description
# -------------------------------------------------
# !			Error lines
# <			Low score alignments
# =			Pairs with more than one best score
# ~			Read pair mapped on the same strand
# ?			Segment length out of reasonable range
# @			Mapped positions overlapped	
#

READ_ERROR = '!'
READ_LOW_SCORE = '<'
READ_MULTI_BEST = '='
READ_WRONG_CHROM = '~'
READ_WRONG_SIZE = '?'
READ_OVERLAPPED = '@'

# sort bam file by qname 
# use pysam sort interface

def SortSam(inBam, outBam):
	pysam.sort("-n", inBam, outBam)

def ScoreByCigar(cigar):
	totalCount = 0
	matchCount = 0
	for tag in cigar :
		if(tag[0] == 0):
			matchCount += tag[1]
		totalCount += tag[1]
	if(totalCount == 0):
		return -1
	return int(round(100.0 * matchCount / totalCount))

def ReadPairScoreByCigar(cigar1, cigar2):
	return ScoreByCigar(cigar1) + ScoreByCigar(cigar2)

def ReadPairLen(read1, read2):
	if(not (read1 and read2)):
		return -1
	return (read2.pos + read2.query_length - read1.pos)

# Unique Single reads
#
# dictSingle structure:
# dictSingle = {"key1": [read1, read2, ...], "key2": [read1, read2, ...]}
#
# Unique Strategy:
#
# 1. For reads with the same qname, keep the alignment with the highest score
# 

def UniqueSingleReads(dictSingle, strictMode, outFile, logFile):
	count = 0

	if(len(dictSingle) <= 0):
		return 0

	for reads in dictSingle.itervalues():
		bestScore = -1
		bestReads = []
		for read in reads:
			score = ScoreByCigar(read.cigar)
			if(score > bestScore):
				bestScore = score
				for bestRead in bestReads :
					
					# write read to log file

					logMsg = READ_LOW_SCORE + ' ' + str(bestRead)
					logFile.write(logMsg + '\n')
				bestReads = [read]
			elif(score == bestScore):
				bestReads += [read]
			else:
				logMsg = READ_LOW_SCORE + ' ' + str(read)
				logFile.write(logMsg + '\n')
		
		# write best reads

		bestReadCount = len(bestReads)
		if(bestReadCount <= 0):
			return 0
		elif(bestReadCount == 1):
			bestRead = bestReads[0]
		else:
			if(strictMode):
				for bestRead in bestReads :
					logMsg = READ_MULTI_BEST + ' ' + str(bestRead)
					logFile.write(logMsg + '\n')
				return 0
			else:
				bestRead = bestReads[0]

		if(bestRead):
			outFile.write(bestRead)
			count += 1

	return count

# Unique paired reads:
#
# dictPaired Structure:
# dictPaired = {'key1': [read1, read2], 'key2':[read1, read2]}
#
# Unique Strategy:
#
# 1. Keep alignment pairs with the highest score.
#
# 2. Reads are properly mapped (on the same chromsome & on the opposite strand)
#
# 3. Reads length should be in a proper range
#

def UniquePairedReads(dictPaired, strictMode, bamFile, outFile, logFile):
	bestScore = -1
	bestPairs = []
	pairCount = 0

	if(len(dictPaired) <= 0):
		return 0

	for key, pair in dictPaired.iteritems():

		# pairs are properly mapped

		if(not pair[0] and not pair[1]):
			logMsg = READ_ERROR + ' ' + key 
			logFile.write(logMsg + '\n')
			continue
		elif(not pair[0] and pair[1]):
			logMsg = READ_ERROR + ' ' + str(pair[1])
			logFile.write(logMsg + '\n')
			continue
		elif(pair[0] and not pair[1]):
			logMsg = READ_ERROR + ' ' + str(pair[0])
			logFile.write(logMsg + '\n')
			continue

		# paired reads on the same chromsome

		if(bamFile.getrname(pair[0].rname) != bamFile.getrname(pair[1].rname)):
			logMsg = READ_WRONG_CHROM + ' ' + str(pair[0])
			logFile.write(logMsg + '\n')
			logMsg = READ_WRONG_CHROM + ' ' + str(pair[1])
			logFile.write(logMsg + '\n')
			continue

		# pair size

		pairLen = ReadPairLen(pair[0], pair[1])
		if(pairLen < lowerLimit or pairLen > upperLimit):
			logMsg = READ_WRONG_SIZE + ' ' + str(pair[0])
			logFile.write(logMsg + '\n')
			logMsg = READ_WRONG_SIZE + ' ' + str(pair[1])
			logFile.write(logMsg + '\n')
			continue

		# pair score 

		score = ReadPairScoreByCigar(pair[0].cigar, pair[1].cigar)
		if(score > bestScore):
			bestScore = score
			for bestPair in bestPairs:
				logMsg = READ_LOW_SCORE + ' ' + str(bestPair[0])
				logFile.write(logMsg + '\n')
				logMsg = READ_LOW_SCORE + ' ' + str(bestPair[1])
				logFile.write(logMsg + '\n')
			bestPairs = [pair]
		elif(score == bestScore):
			bestPairs += [pair]
		else:
			logMsg = READ_LOW_SCORE + ' ' + str(pair[0])
			logFile.write(logMsg + '\n')
			logMsg = READ_LOW_SCORE + ' ' + str(pair[1])
			logFile.write(logMsg + '\n')

	# write best pair

	bestPairCount = len(bestPairs)
	if(bestPairCount <= 0):
		return 0
	elif(bestPairCount == 1):
		bestPair = bestPairs[0]
	else:
		if(strictMode):
			for bestPair in bestPairs:
				logMsg = READ_MULTI_BEST + ' ' + str(bestPair[0])
				logFile.write(logMsg + '\n')
				logMsg = READ_MULTI_BEST + ' ' + str(bestPair[1])
				logFile.write(logMsg + '\n')
			return 0
		else:
			bestPair = bestPairs[0]

	if(bestPair):
		outFile.write(bestPair[0])
		outFile.write(bestPair[1])
		pairCount = 2

	return pairCount

def main():
	# parse the command line options
	usage = 'usage: %prog [options] input.bam -o output.bam'
	parser = OptionParser(usage=usage, version='%prog version 0.1.0')
	parser.add_option('-o', '--output-file', dest='outputfile',
						help='write the result to output file')
	parser.add_option('-s', '--sort', 
						action="store_true", dest="sort", default=False,
						help='sort the input SAM file before further processing')
	parser.add_option('-k', '--reg-key', dest='regkey',
						help='regexp to extract the key part of the read qname')
	parser.add_option('-n', '--reg-num', dest='regnum',
						help='regexp to extract the number part of the read qname')
	parser.add_option('-u', '--upper-limit', dest='upperlimit',
						help='the upper limit of the paired reads length')
	parser.add_option('-l', '--lower-limit', dest='lowerlimit', 
						help='the lower limit of the paired reads length')
	parser.add_option('-t', '--strict-mode', 
						action="store_true", dest="strictmode", default=False,
						help='strict mode')
	parser.add_option('-q', '--quiet-mode', 
						action="store_true", dest="quiet", default=False,
						help='supress progress output')

	(options, args) = parser.parse_args()
	if(len(args) != 1):
		parser.print_help()
		sys.exit(0)
	inputBamFileName = args[0]
	baseFileName = os.path.splitext(os.path.basename(inputBamFileName))[0]
	bamFileName = inputBamFileName
	rmTemp = False

	if(options.sort):
		print('[*] Sorting by QNAME...')
		bamFileName = 'sorted.' + baseFileName
		SortSam(inputBamFileName, bamFileName)
		bamFileName += '.bam'
		rmTemp = True

	strictMode = False
	if(options.strictmode):
		strictMode = True

	# Prepare the qname regexp

	global readKeyRe
	global readNumRe
	if(options.regkey):
		readKeyRe = re.compile(options.regkey)
	if(options.regnum):
		readNumRe = re.compile(options.regnum)

	global upperLimit
	global lowerLimit
	if(options.upperlimit):
		upperLimit = int(options.upperlimit)
	if(options.lowerlimit):
		lowerLimit = int(options.lowerlimit)

    # Load the input bam file

	print('[*] Initializing...')

	if(not os.path.exists(bamFileName)):
		print('error: Failed to open file "', bamFileName, '"')
		sys.exit(-1)
	bamFile = pysam.AlignmentFile(bamFileName, "rb")
	
	# prepare the output file

	outputBamFileName = 'unique.' + baseFileName + '.bam'
	if(options.outputfile):
		outputBamFileName = options.outputfile

	outFile = pysam.AlignmentFile(outputBamFileName, 'wb', template = bamFile)

	logFileName = baseFileName + '.rmdup.log'
	try:
		logFile = open(logFileName, 'w')
	except IOError:
		print('error: Create log file failed!')
		sys.exit(-1)

	# process file

	readCount = 0
	keepedCount = 0
	print('[*] Processing...')
	dictPaired = {}
	dictSingle = {}
	currentGroupKey = ''
	for read in bamFile.fetch(until_eof = True):
		try:
			groupKey = ''.join(readKeyRe.findall(read.qname)[0])
		except:
			groupKey = read.qname

		if(groupKey != currentGroupKey):
			currentGroupKey = groupKey

			# handle and write results

			keepedCount += UniquePairedReads(dictPaired, strictMode, bamFile, outFile, logFile)
			keepedCount += UniqueSingleReads(dictSingle, strictMode, outFile, logFile)
			dictPaired.clear()
			dictSingle.clear()
		
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

			if(readKey in dictPaired):
				if(dictPaired[readKey][readType]):
					logMsg = READ_OVERLAPPED + ' ' + str(read)
					logFile.write(logMsg + '\n')
				else:
					dictPaired[readKey][readType] = read
			else:
				dictPaired[readKey] = [None, None]
				dictPaired[readKey][readType] = read
		else:

			# single read

			try:
				readKey = groupKey + ':' + readNumRe.findall(read.qname)[0]
			except:
				readKey = groupKey

			if(readKey in dictSingle):
				dictSingle[readKey] += [read]
			else:
				dictSingle[readKey] = [read]

	
		# progress information 

		readCount = readCount + 1
		if not options.quiet :
			sys.stdout.write('\r    read #%ld' % (readCount))
			sys.stdout.flush()

	# write the cached alignments

	keepedCount += UniquePairedReads(dictPaired, strictMode, bamFile, outFile, logFile)
	keepedCount += UniqueSingleReads(dictSingle, strictMode, outFile, logFile)

	if(readCount == 0):
		keepedPercent = 0.0
	else:
		keepedPercent = keepedCount * 1.0 / readCount
	print('\n    %ld (%.2f%%) alignments written.' % (keepedCount, keepedPercent * 100))
	
	# Clear resources

	bamFile.close()
	outFile.close()
	logFile.close()

	if(rmTemp):
		os.unlink(bamFileName)

	print('[*] Complete')

if __name__ == '__main__':
	main()
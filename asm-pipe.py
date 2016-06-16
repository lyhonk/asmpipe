#! /usr/bin/env  python
"""
Usage:
  asmpipe.py [-rsk] CONFIGFILE BAMFILE
  asmpipe.py --version

Arguments:
  CONFIGFILE     asm configfile
  BAMFILE        input.bam

Options:
  -r --rmdup       remove duplicate read, usually chose options -rsk together
  -s --sort        sort the input bamfile
  -k --keep        keep the sorted bamfile
  -h --help        show the manual
  -v --version     show the version
"""

try:
	from docopt import docopt
except ImportError:
	exit('This script requires package `docopt` \n'
	'install this package first: sudo pip install docopt\n')

import ConfigParser
import string, os, sys, getpass, logging, shutil
sys.path.append("module")
from rmdup import *

def generateTrakHub( fileBaseName, trackhubPath, genomeType, tagMethBw, tagUmethBw, tagAsmBb):
	trackGenomePath = trackhubPath + "/" + genomeType
	os.makedirs(trackGenomePath)
	trakDataPath =  trackGenomePath + "/data"
	os.makedirs(trakDataPath) 
	shutil.move(tagMethBw, trakDataPath)
	shutil.move(tagUmethBw, trakDataPath)
	shutil.move(tagAsmBb, trakDataPath)
	## hub. txt ##
	hubFile = trackhubPath + "hub.txt"
	emailAddr = getpass.getuser() + "@" + os.uname()[1]
	f=open("tpl/hub.txt", 'r').readlines()
	for i in range( len(f) - 1 ):
		f[i]=f[i].replace("$email$",emailAddr)
	tmpFile = open(hubFile, 'w')
	tmpFile.writelines(f)
	tmpFile.close()
	## genomes.txt ##
	genomesFile = trackhubPath + "/genomes.txt"
	f = open(genomesFile,"w")
	f.writelines("genome " + genomeType + "\n")
	f.writelines("trackDb " + genomeType + "/trackDb.txt\n")
	f.close()
	## hub.html ##
	hubHtml = trackhubPath + "/hub.html"
	f = open("tpl/hub.html", 'r').readlines()
	for i in range( len(f) - 1 ):
		f[i]=f[i].replace("$genomeType$",genomeType)
		f[i]=f[i].replace("$description$", "description")
	tmpFile = open(hubHtml, 'w')
	tmpFile.writelines(f)
	tmpFile.close()
	## trackDb.txt ##
	trackDbFile = trackhubPath + "/trackDb.txt"
	f = open( "tpl/trackDb.txt", "r").readlines()
	for i in range( len(f) - 1 ):
		f[i] = f[i].replace( "$ASMtrackName$", fileBaseName + ".asm")
		f[i] = f[i].replace( "$asmURL$", "./data/" + tagAsmBb)
		f[i] = f[i].replace( "$asmLabel$", fileBaseName + ".asm")
		f[i] = f[i].replace( "$methName$", fileBaseName + ".meth")
		f[i] = f[i].replace( "$methLabel$", fileBaseName + ".meth")
		f[i] = f[i].replace( "$methURL$", "./data/" + tagMethBw)
		f[i] = f[i].replace( "$umethName$", fileBaseName +".umeth" )
		f[i] = f[i].replace( "$umethLabel$", fileBaseName + ".umeth" )
		f[i] = f[i].replace( "$umethURL$", "./data/" + tagUmethBw)
	tmpFile = open(trackDbFile, 'w')
	tmpFile.writelines(f)
	tmpFile.close()


if __name__ == '__main__':
	arguments = docopt(__doc__, version='Version 0.0.1')
	configPath = arguments['CONFIGFILE']
	bamPath = arguments['BAMFILE']
	if( not( os.path.exists(configPath) & os.path.exists(bamPath) ) ):
		print ("error: %s or %s doest not exist."%(configPath, bamPath)) 
		sys.exit(-1)

	## tool cmds ##
	scriptPath = os.getcwd()
	chromsizePath = scriptPath + "/chrom.sizes"
	wigtobigwig = scriptPath + "/utils/wigToBigWig"
	bedtobigbed = scriptPath + "/utils/bedToBigBed"

	## file path and names ##
	inputBamBasename = os.path.basename(bamPath)
	inputBamBasename = os.path.splitext(inputBamBasename)[0]
	logFile = inputBamBasename + ".pipeline.log"
	progressFile = "." + inputBamBasename + ".progress"

	## read config file ##
	cf = ConfigParser.ConfigParser() 
	cf.read(configPath) 
	refFaFile = cf.get("conf","ref_Fa")
	genomeType = cf.get("conf","genomeType")
	if( not( os.path.exists(refFaFile) ) ):
		print(" Invalid reference file ref_fa=%s in %s "%( refFaFile, configPath) )
		sys.exit( -1 )
	if ( genomeType not in ["mm9", "mm10", "hg18", "hg19"] ):
		print(" Error: unknown genome_type - %s, only mm9 | mm10 | hg18 | hg19 are surported"%(genomeType))
		sys.exit( -1 )

	##  Initialize ##
	# set logger #
	logging.basicConfig(level=logging.DEBUG,
		format='%(asctime)s %(message)s',
		datefmt='%Y-%m-%d %H:%M:%S',
		filename=logFile,
		filemode='w')
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	formatter = logging.Formatter('%(message)s')
	console.setFormatter(formatter)
	logging.getLogger('').addHandler(console)
	logger = logging

	logger.info("[*] Preparing output files & path...")

	# prepare progress file #
	if (os.path.exists(progressFile)):
		progLog = ConfigParser.ConfigParser()
		progLog.read(progressFile)
	else:
		os.mknod(progressFile) 
		progLog = ConfigParser.ConfigParser()
		progLog.add_section("progress")
		progLog.write(open(progressFile, "w"))

	# output path #
	outputPath = inputBamBasename + ".tagmeth"	
	if ( not ( os.path.exists(outputPath) ) ) : os.mkdir( outputPath)

	trackhubPath = inputBamBasename + ".trackhub"
	if ( not ( os.path.exists( trackhubPath) ) ) : os.mkdir( trackhubPath)

	logger.info("    Everything has been prepared !")

	## setp1 : clean up reads ##
	if(arguments['--rmdup']):
		logger.info("[*] Cleaning up reads...")

		try:
			progLog.get("progress","asmCleanRead")
			logger.warning("\n[***] Pay attation !  We found rmdup progress done before !\n")
			pass

		except ConfigParser.NoOptionError:
			keyRe = cf.get("conf", "keyRe")
			numRe = cf.get("conf", "numRe")
			readKeyRe = re.compile(keyRe)
			readNumRe = re.compile(numRe)
			rmdup = rmdup(bamPath, sortBam = arguments['--sort'], keepSorted = arguments['--keep'])
			rmdupPath = inputBamBasename + ".rmdup"
			if ( not ( os.path.exists( rmdupPath ) ) ) : os.mkdir( rmdupPath )
			shutil.move( inputBamBasename + ".rmdup.log", rmdupPath) 
			shutil.move("unique." + inputBamBasename +".bam", rmdupPath)
			if ( arguments['--keep'] ): shutil.move( "sorted." + inputBamBasename + ".bam", rmdupPath)
			progLog.set("progress","asmCleanRead","OK")
			progLog.write(open(progressFile, "w"))
			logger.info(rmdup)
			logger.info("[*] " + bamPath + " has removed duplicated reads")

	## setp2 : generate tag meth csv file##



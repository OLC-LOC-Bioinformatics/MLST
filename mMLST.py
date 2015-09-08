#! /usr/env/python
__author__ = 'mikeknowles, akoziol'
""" Includes threading found in examples:
http://www.troyfawkes.com/learn-python-multithreading-queues-basics/
http://www.ibm.com/developerworks/aix/library/au-threadingpython/
https://docs.python.org/2/library/threading.html
Revised with speed improvements
"""
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from threading import Thread
from Queue import Queue
from collections import defaultdict
from cStringIO import StringIO
from glob import glob
import subprocess, os, time, sys, shlex, re, threading, json, errno, operator
from argparse import ArgumentParser

parser = ArgumentParser(description='Performs blast analyses to determine presence of alleles in a genome query, and '
                                    'types genome based on typing profile. Adds novel alleles and profiles to the '
                                    'appropriate files. '
                                    'Example command: '
                                    '-p /home/blais/PycharmProjects/MLST  '
                                    '-s /home/blais/PycharmProjects/MLST/sequences '
                                    '-O /home/blais/PycharmProjects/MLST/Organism '
                                    '-o Vibrio '
                                    '-S MLST')
parser.add_argument('-p', '--path', required=True,
                    # default='/home/blais/PycharmProjects/pythonGeneSeekr/',
                    help='Specify path for custom folder locations. If you don\'t supply additional paths'
                         'e.g. sequencePath, allelePath, or organismPath, then the program will look for '
                         'MLST files in .../path/Organism, and the query sequences in ../path/sequences')
parser.add_argument('-c', '--cutoff', required=False, default=100,
                    help='The percent identity cutoff value for BLAST matches. Default is 100%)')
parser.add_argument('-s', '--sequencePath', required=False,
                    default='/home/blais/PycharmProjects/MLST/sequences',
                    help='The location of the query sequence files')
parser.add_argument('-a', '--alleleProfilePath', required=False,
                    # default='/home/blais/PycharmProjects/pythonGeneSeekr/Organism/Salmonella/cgMLST',
                    help='The path of the folder containing the two folders containing '
                         'the allele files, and the profile file e.g. /path/to/folder/Organism/Salmonella/cgMLST')
parser.add_argument('-O', '--organismPath', required=False,
                    help='The path of the folder containing the organism folders e.g. /path/to/folder/Organism')
parser.add_argument('-o', '--organism', required=False,
                    help='The name of the organism you wish to type. Must match the folder name containing the schemes'
                         'e.g. Salmonella')
parser.add_argument('-S', '--scheme', required=False,
                    help='The scheme you wish to use. Must match the folder name containing the scheme e.g. cgMLST.'
                         'Furthermore, this folder must contain two folders: "alleles" and "profile". The alleles '
                         'folder contains the allele files in .fasta format, and the profile folder contains '
                         'the profile in .txt format.')
parser.add_argument('-u', '--updateProfileFalse', required=False, default=True,
                    help='By default, the program automatically creates new sequence profiles and appends these '
                         'profiles to the profile file. If, instead, you wish to wish to see the closest match of a '
                         'query genome to known reference profiles, set this to False.')

# Get the arguments into a list
args = vars(parser.parse_args())

# Define variables from the arguments - there may be a more streamlined way to do this
# Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
path = os.path.join(args['path'], "")
cutoff = float(args['cutoff'])/100
if args['sequencePath']:
    sequencePath = os.path.join(args['sequencePath'], "")
else:
    sequencePath = ""
if args['alleleProfilePath']:
    allelePath = os.path.join(args['alleleProfilePath'], "")
else:
    allelePath = ""
if args['organismPath']:
    organismPath = os.path.join(args['organismPath'], "")
else:
    organismPath = ""
scheme = args['scheme']
organism = args['organism']
updateProfile = args['updateProfileFalse']

# Empty the updateProfile if it is not True - will be used in if checks later
if updateProfile != True:
    updateProfile = ""


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
        # os.chmod(inPath, 0777)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    return defaultdict(make_dict)

# Initialise the count used in the dotter function
count = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # Use a global variable
    global count
    if count <= 80:
        sys.stdout.write('.')
        count += 1
    else:
        sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        count = 1


def makeblastdb(dqueue):
    """Makes blast database files from targets as necessary"""
    while True:  # while daemon
        fastapath = dqueue.get()  # grabs fastapath from dqueue
        # remove the path and the file extension for easier future globbing
        db = fastapath.split(".")[0]
        nhr = "%s.nhr" % db  # add nhr for searching
        # print nhr
        FNULL = open(os.devnull, 'w')  # define /dev/null
        if not os.path.isfile(str(nhr)):  # if check for already existing dbs
            subprocess.Popen(shlex.split("makeblastdb -in %s -dbtype nucl -out %s" % (fastapath, db)), stdout=FNULL, stderr=FNULL)
            # make blastdb
            dotter()
        dqueue.task_done()  # signals to dqueue job is done
        sys.exit()  # necessary

# Declare queues, list, and dictionaries
dqueue = Queue()
blastqueue = Queue()
testqueue = Queue()
plusqueue = Queue()
plusdict = defaultdict(make_dict)
profileData = defaultdict(make_dict)
MLSTseqType = defaultdict(make_dict)
bestDict = defaultdict(make_dict)
resultProfile = defaultdict(make_dict)
genedict = defaultdict(list)
blastpath = []
threadlock = threading.Lock()


def makedbthreads(fastas):
    """Setup and create threads for class"""
    # Create and start threads for each fasta file in the list
    for i in range(len(fastas)):
        # Send the threads to makeblastdb
        threads = Thread(target=makeblastdb, args=(dqueue,))
        # Set the daemon to true - something to do with thread management
        threads.setDaemon(True)
        # Start the threading
        threads.start()
    for fasta in fastas:
        # Add the fasta file to the queue
        dqueue.put(fasta)
    dqueue.join()  # wait on the dqueue until everything has been processed


def xmlout(fasta, genome):
    """Parses variables from supplied tuples? dictionaries?"""
    global path
    # Extract the variables from the passed variables
    gene = fasta.split('/')[-1]  # split file from path, could use os.path.split
    genename = gene.split('.')[0]
    genomename = os.path.basename(genome).split(".")[0]
    # Create the tmp folder (if necessary)
    make_path("%stmp" % path)
    out = "%stmp/%s.%s.xml" % (path, genomename, genename)  # Changed from dictionary to tuple
    # Return the parsed variables
    return path, gene, genename, genomename, out, fasta


def alleleSplitter(alleleNames):
    # Multiple try-excepts. Maybe overly complicated, but I couldn't get it work any other way
    # This (hopefully) accounts for all the possible naming schemes for the alleles
    try:  # no split - just allele numbers e.g. >12
        match = re.search(r"(>\d+)", alleleNames)
        alleleNumber = match.group().split(">")[1]
        allelePreNumber = ""
    except (IndexError, AttributeError):
        try:  # split on "_" e.g. >AROC_12
            # alleleNumber is the number of the allele(!). It should be different for each allele
            alleleNumber = alleleNames.split("_")[1]
            # allelePreNumber is anything before the allele number. It should be the same for each allele
            allelePreNumber = alleleNames.split("_")[0]
        except IndexError:
            try:  # split on "-" e.g. >AROC-12
                alleleNumber = alleleNames.split("-")[1]
                allelePreNumber = alleleNames.split("-")[0]
            except IndexError:
                try:  # split on " " e.g. >AROC 12
                    alleleNumber = alleleNames.split(" ")[1]
                    allelePreNumber = alleleNames.split(" ")[0]
                except IndexError:
                    try:  # split on change from letters to numbers e.g. >AROC12
                        match = re.search(r"(>[A-z/a-z]+)(\d+)", alleleNames)
                        alleleNumber = match.groups()[1]
                        allelePreNumber = match.groups()[0]
                    except (IndexError, AttributeError):
                        alleleNumber = alleleNames
                        allelePreNumber = alleleNames
    # Return the variables
    return int(alleleNumber), allelePreNumber


def blastparse(blast_handle, genome, gene, cutoff, genepath):
    """Parses BLAST results, and populates a dictionary with the results"""
    global plusdict
    global profileData
    snpDict = {}
    dataDict = {}
    records = NCBIXML.parse(blast_handle)   # Open record from memory-mapped file
    dotter()
    # Split the extension from the genome name
    genomeName = os.path.basename(genome).split('.')[0]
    numhsp = sum(line.count('<Hsp>') for line in iter(blast_handle.readline, ""))
    if numhsp >= 1:
        # Since we scanned through result_handle looking for HSPs, the position of the read/write pointer
        # within the file is at the end. To reset it back to the beginning, .seek(0) is used
        blast_handle.seek(0)
        # There was an issue with errors referring to end of record - I think this had something to do with the program
        # trying to read .xml files that improperly formatted due to the program crashing mid-file creation. I don't
        # know if this actually does anything, or whether there just haven't been any issues since I made this change
        # if records:
        for record in records:  # This process is just to retrieve HSPs from xml files
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    # Calculate the percent identity
                    percentIdentity = "%.2f" % float(float(hsp.identities) / float(alignment.length) * 100)
                    allele = str(alignment.title.split(" ")[-1])
                    # If the results are greater than the cutoff value, add them to the dictionary
                    if hsp.identities >= alignment.length * cutoff:
                        # Clears out any "N" values in the dictionary
                        if "N" in plusdict[genomeName][gene]:
                           plusdict[genomeName][gene].clear()
                        plusdict[genomeName][gene][allele] = percentIdentity
                        # As the blast results files are not sorted by percent identity, and, at least for rMLST
                        # genes, not all genes are the same length, a match can occur after lots of close matches
                        snpDict.clear()
                    # Process results if the percent identity is 98-99.9% and there is no 100% match in the dictionary
                    # elif not gene in plusdict[genomeName] and hsp.identities >= alignment.length * cutoff * 0.98:
                    elif not gene in plusdict[genomeName] and not snpDict and hsp.identities >= alignment.length * cutoff * 0.98:
                        # Puts the HSP in the correct order -  hits to the negative strand will be
                        # reversed compared to what we're looking for
                        if hsp.sbjct_start < hsp.sbjct_end:
                            end = hsp.sbjct_end
                        else:
                            end = hsp.sbjct_start
                        # Screen out hits that are shorter than the targets
                        # Keeping it this format even though this if statement could be re-written more efficiently
                        if end < alignment.length:
                            pass
                        else:
                            # Add the details of the mismatching allele to two dictionaries to be processed below
                            snpDict[genepath] = hsp.query
                            dataDict[genomeName] = gene
                    # If the percent identity is below the 98% cutoff threshold or is the hsp is too short
                    elif not gene in plusdict[genomeName] and not snpDict:
                        plusdict[genomeName][gene]["N"] = 0
    # If there are no records, populate the dictionary with "N" and a 0% identity
    else:
        plusdict[genomeName][gene]["N"] = 0
    # Add matches that are 98.0 < match identity < 99.9 and the same length of the target to the allele file
    if snpDict:
        # Initialise some variables
        alleleNames = []  # The allele names already in the allele file
        allelePreNumber = ""  # The text before the allele number e.g. >AROC
        alleleNumber = 0  # The last allele in the file e.g. 72
        # Go through the allele files in snpDict
        for gPath in snpDict:
            # Open the allele file
            with open(gPath) as geneFile:
                for line in geneFile:
                    # Only interested in the header for each allele
                    if ">" in line:
                        # Append all, but only use the last - used to be a string instead of a list
                        alleleNames.append(line)
            # Find the allele number and the text before the number for different formats
            alleleNumber, allelePreNumber = alleleSplitter(alleleNames[-1])
            # I wanted to keep the new allele numbers distinct from the old ones, so they will start at 1000000
            if alleleNumber < 1000000:
                newAlleleNumber = 1000000
            # As this will continuously update the allele database, I need to check if I've already added new alleles
            else:
                newAlleleNumber = alleleNumber + 1
            # Initialise newAllele - formatted updated allele number e.g. >AROC1000000
            newAllele = ""
            # Accounts for no pre-number text being present (e.g. >12)
            if not allelePreNumber:
                # Create a sequence record using BioPython
                fasta = SeqRecord(Seq(snpDict[gPath], "fasta"),  # snpDict[gPath] is the hsp sequence
                                  description="",  # if this is not added, then some junk is written to the new header
                                  id=str(newAlleleNumber))  # keep the format of just the allele number e.g. >1000000
            # If there is pre-number text, do the same thing, but use the allele pre-number as well
            else:
                allelePreNumber = allelePreNumber.split(">")[1]
                newAllele = "%s-%s" % (allelePreNumber, newAlleleNumber)
                fasta = SeqRecord(Seq(snpDict[gPath], "fasta"),
                                  description="",
                                  id=newAllele)
            # Open the allele file to append the new sequence record
            with open(gPath, "a") as updatedFasta:
                # Use the SeqIO module to properly format the new sequence record
                SeqIO.write(fasta, updatedFasta, "fasta")
            # Perform some cleanup - need to remove all database and results files associated with the allele file
            # prior to its update. Otherwise, the new changes to the file will not be reflected on subsequent iterations
            # of this genome/allele file comparison
            # The path and file name of the allele file without an extension
            baseName = os.path.splitext(gPath)[0]  # Despite the fact that I don't use os.path.basename to generate this
            # Remake the database files
            updatedb = [gPath]
            makedbthreads(updatedb)

            # # Try to remove all the appropriate files
            # try:
            #     os.remove("%s.nhr" % baseName)
            #     os.remove("%s.nin" % baseName)
            #     os.remove("%s.nsq" % baseName)
            #
            # except OSError:
            #     raise
            # Now use this new allele in populating plusdict
            for updatedGenome in dataDict:
                # The percent identity has to be 100% - this allele matches itself
                plusdict[updatedGenome][dataDict[updatedGenome]][newAllele] = 100.00
    # Should this have been closed earlier? I don't know
    blast_handle.close()


class runblast(threading.Thread):
    """Runs the multi-threaded blast analysis"""
    def __init__(self, blastqueue):
        self.blastqueue = blastqueue
        threading.Thread.__init__(self)

    def run(self):
        while True:
            global blastpath, plusdict  # global varibles, might be a better way to pipe information
            genome, fasta, blastexist, cutoff = self.blastqueue.get()  # retrieve variables from queue
            path, gene, genename, genomename, out, genepath = xmlout(fasta, genome)  # retrieve from string splitter
            #Precaution
            threadlock.acquire()
            # Add the appropriate variables to blast path
            blastpath.append((out, path[-1], gene, genename,))  # tuple-list
            threadlock.release()
            # Print a dot for each gene, genome combination processed
            dotter()
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db,
            # a mild evalue of 0.1, and XML formatted output
            # Removed perc_identity=percentIdentity from the call, as this allows more flexibility for parsing files
            # with different cutoff values - if I want to re-analyse a search with a lower cutoff, I won't have to
            # re-perform the BLAST search each time
            db = fasta.split(".")[0]
            blastn = NcbiblastnCommandline(query=genome, db=db, evalue=0.1, outfmt=5)
            # Note that there is no output file specified -  the search results are currently stored in stdout
            stdout, stderr = blastn()
            # Search stdout for matches - if the term Hsp appears (the .find function will NOT
            # return -1), a match has been found, and stdout is written to file
            if stdout.find('Hsp') != -1:
                blast_handle = StringIO(stdout)  # Convert string to IO object for use in SearchIO using StringIO
                blastparse(blast_handle, genome, genename, cutoff, genepath)  # parse the data already in memory
            # If there are no hsps, then populate the dictionary with the negative results
            else:
                plusdict[genomename][genename]["N"] = 0
            # Close the thread
            self.blastqueue.task_done()


def blastnthreads(fastas, genomes, cutoff):
    """Setup and create  threads for blastn and xml path"""
    blastexist = {}
    # Create threads for each gene, genome combination
    for genome in genomes:
        for fasta in fastas:
            # Add the appropriate variables to the threads
            blastqueue.put((genome, fasta, blastexist, cutoff))
        blastqueue.join()


def blastDatabaseClearer(genePath):
    """Due to the nature of the program updating allele files, it's not desirable to use previously generated databases.
    Additionally, with the use of these files by multiple programs, there is an issue. This script makes database files
    as follows: aroC.fasta becomes aroC.nhr, etc. The current SPAdes assembly pipeline would take that same .fasta file
    and create aroC.fasta.nhr. Deleting database files prevents issues with glob including database files."""
    # Get all the .nhr, .nin, .nsq files
    databaseList = glob("%s/*.n*" % genePath)
    # And delete them
    for allele in databaseList:
        os.remove(allele)


def organismChooser(path, alleleProfilePath, organismPath, schemeName, organismName):
    """Allows the user to choose which organism to be used in the analyses"""
    # Initialise a count variable to be used in extracting the desired entry from a list of organisms
    orgcount = 0
    schemecount = 0
    alleles = []
    profile = []
    # If the path of the folder containing the allele and profile subfolders is provided
    if alleleProfilePath:
        # Set the genePath variable for use in blastDatabaseClearer
        genePath = "%salleles" % alleleProfilePath
        # Remove and previously created blast database files
        blastDatabaseClearer(genePath)
        # Create lists of the alleles, and the profile
        alleles = glob("%salleles/*.*fa*" % alleleProfilePath)
         # If the profile has previously been processed in this script, use the .json file
        profile = glob("%sprofile/*.json" % alleleProfilePath)
        # Otherwise use the .txt file
        if not profile:
            profile = glob("%sprofile/*.txt" % alleleProfilePath)
    else:
        # Get a list of the organisms in the (default) Organism subfolder
        if not organismPath and not organismName:
            orgList = glob("%sOrganism/*" % path)
            # Iterate through the sorted list
            for folder in sorted(orgList):
                # Ensure that folder is, in actuality, a folder
                if os.path.isdir(folder):
                    # Print out the folder names and the count
                    print "[%s]: %s" % (orgcount, os.path.split(folder)[1])
                    orgcount += 1
            # Get the user input - the number entered corresponds to the list index
            response = input("Please select an organism: ")
            # Get the organism path into a variable
            organism = sorted(orgList)[int(response)]
            organismName = os.path.split(organism)[1]
        elif organismPath and not organismName:
            orgList = glob("%s*" % organismPath)
            # Iterate through the sorted list
            for folder in sorted(orgList):
                # Ensure that folder is, in actuality, a folder
                if os.path.isdir(folder):
                    # Print out the folder names and the count
                    print "[%s]: %s" % (orgcount, os.path.split(folder)[1])
                    orgcount += 1
            # Get the user input - the number entered corresponds to the list index
            response = input("Please select an organism: ")
            # Get the organism path into a variable
            organism = sorted(orgList)[int(response)]
            organismName = os.path.split(organism)[1]
        # Get the schemes into a list
        if not schemeName:
            if not organismPath:
                schemeList = glob("%sOrganism/%s/*" % (path, organismName))
            else:
                schemeList = glob("%s/%s/*" % (organismPath, organismName))
            # Iterate through the sorted list
            for folder in sorted(schemeList):
                # Ensure that folder is, in actuality, a folder
                if os.path.isdir(folder):
                    # Print out the folder names and the count
                    print "[%s]: %s" % (schemecount, os.path.split(folder)[1])
                    schemecount += 1
            # Same as above
            schemeResponse = input("Please select a typing scheme:")
            scheme = sorted(schemeList)[int(schemeResponse)]
            schemeName = os.path.split(scheme)[1]
            # Set the variables as above
            genePath = "%s/alleles" % scheme
            blastDatabaseClearer(genePath)
            alleles = glob("%s/alleles/*.*fa*" % scheme)
            profile = glob("%s/profile/*.json" % scheme)
            if not profile:
                profile = glob("%s/profile/*.txt" % scheme)
        # If the name of the typing scheme is provided
        else:
            # If the organism path is not provided
            if not organismPath:
                # Default to using "Organism" in the path
                scheme = "%sOrganism/%s/%s" % (path, organismName, schemeName)
            else:
                # Otherwise set scheme as follows:
                scheme = "%s%s/%s" % (organismPath, organismName, schemeName)
            # Put the query and quality genes into lists
            genePath = "%s/alleles" % scheme
            blastDatabaseClearer(genePath)
            alleles = glob("%s/alleles/*.*fa*" % scheme)
            profile = glob("%s/profile/*.json" % scheme)
            if not profile:
                profile = glob("%s/profile/*.txt" % scheme)
    # Get the path for the original .txt profile file
    profileTxt = "%s.txt" % os.path.splitext(profile[0])[0]
    return alleles, profile[0], organismName, schemeName, profileTxt


def profilR(profileFile):
    """Creates a dictionary from the profile scheme"""
    # Initialise the dictionary
    global profileData
    lastEntry = ""
    geneList = []
    # The gene names are present in the first line of the profile file
    # Note: if the profile files are ever updated, then the clonal complex column must be removed
    # Make the json filename from profileFile - it might already be .json, but this doesn't take long to do
    JSONProfile = "%s.json" % os.path.splitext(profileFile)[0]
    # If this scheme has previously been used, then the profileData dictionary is written to disk for increased speed.
    # Parsing a json file was approximately 10 times faster than parsing the original tab-delimited file
    # Get the MLST profiles for each sequence type
    with open(profileFile) as profile:
        # Files have to be in tab-delimited format
        header = profile.readline().rstrip().split("\t")
        # As certain typing schemes are not in any discernible order, using a naturally ordered list instead of a
        # dictionary to store the names is a good idea
        for gene in header:
            # The first column must have "ST" in the header
            if not "ST" in gene:
                dotter()
                geneList.append(gene)
        # Don't do this if the .json profile has previously been created
        if not os.path.isfile(JSONProfile):
            for line in profile:
                # Grab the name of the last profile
                # MLSTcount will used to associate the gene name in header to the allele (e.g. adk 12)
                MLSTcount = 1
                # Don't need to get the header information again
                # if not "ST" in line:
                    # len(header) will be the same length as the data in line
                while MLSTcount < len(header):
                    # Remove newlines and split on tabs
                    data = line.rstrip().split("\t")
                    # Populate profileData with the sequence type, gene name, and the allele number
                    profileData[data[0]][header[MLSTcount]] = data[MLSTcount]
                    # Increment
                    MLSTcount += 1
            # Split the name (if necessary) to just have the profile number
            # Write the json file to disk
            JSONreport = open(JSONProfile, "wb")
            output = json.dumps(profileData, sort_keys=True, indent=4, separators=(',', ': '))
            JSONreport.write(output)
            JSONreport.close()
        else:
            with open(JSONProfile, "rb") as jsonReport:
                # Load the data
                profileData = json.load(jsonReport)
    return profileData, geneList


def sequenceTyper(profileDict, profileFile, geneList, updateProfile):
    """Determines the sequence type of each strain based on comparisons to sequence type profiles"""
    global MLSTseqType
    global bestDict
    # Initialise variables
    header = 0
    alleleCount = 0
    multiAllele = []
    multiPercent = []
    bestMatch = defaultdict(int)
    bestCount = 0
    # Iterate through the genomes
    for genome in plusdict:
        global resultProfile
        # Initialise bestMatch[genome] with an integer - this will eventually be replaced by the number of matches
        bestMatch[genome] = defaultdict(int)
        # For each gene in plusdict[genome]
        for gene in plusdict[genome]:
            # Clear the appropriate count and lists
            alleleCount = 0
            multiAllele = []
            multiPercent = []
            for allele, percentID in plusdict[genome][gene].iteritems():
                # "N" alleles screw up the allele splitter function
                if allele != "N":
                    # Use the alleleSplitter function to get the allele number
                    alleleNumber, allelePreNumber = alleleSplitter(allele)
                    # Append as appropriate - alleleNumber is treated as an integer for proper sorting
                    multiAllele.append(int(alleleNumber))
                    multiPercent.append(percentID)
                # If the allele is "N"
                else:
                    # Append "N" and a percent identity of 0
                    multiAllele.append("N")
                    multiPercent.append(0)
                # Trying to catch cases that where the allele isn't "N", but can't be parsed by alleleSplitter
                if not multiAllele:
                    multiAllele.append("N")
                    multiPercent.append(0)
                # Populate bestDict with genome, gene, alleles - joined with a space (this was written like this because
                # allele is a list generated by the .iteritems() above, and the percent identity
                bestDict[genome][gene][" ".join(str(allele) for allele in sorted(multiAllele))] = multiPercent[0]
            # Find the profile with the most alleles in common with the query genome
            for sequenceType in profileDict:
                # Reset counts to 0
                matchCount = 0
                bestCount = 0
                # The number of genes in the analysis
                header = len(profileData[sequenceType])
                # refallele is the allele number of the sequence type
                refAllele = profileData[sequenceType][gene]
                # If there are multiple allele matches for a gene in the reference profile e.g. 10 692
                if len(refAllele.split(" ")) > 1:
                    # Map the split (on a space) alleles as integers - if they are treated as integers,
                    # the alleles will sort properly
                    intRefAllele = map(int, refAllele.split(" "))
                    # Create a string of the joined, sorted allleles
                    sortedRefAllele = " ".join(str(allele) for allele in sorted(intRefAllele))
                else:
                    # Use the reference allele as the sortedRefAllele
                    sortedRefAllele = refAllele
                for allele, percentID in bestDict[genome][gene].iteritems():
                    # If the allele for the gene in the query genome matches the allele in the reference profile,
                    # add the result to the bestMatch dictionary. Because genes with multiple alleles were sorted the
                    # same way, these strings with multiple alleles will match e.g. 10 692 will never be 692 10
                    if allele == sortedRefAllele:
                        # Increment the number of matches to each profile
                        bestMatch[genome][sequenceType] += 1
        # Get the best number of matches
        # From: https://stackoverflow.com/questions/613183/sort-a-python-dictionary-by-value
        sortedMatches = sorted(bestMatch[genome].items(), key=operator.itemgetter(1), reverse=True)[0][1]
        # If there are fewer matches than the total number of genes in the typing scheme
        if int(sortedMatches) < header:
            # Iterate through the sequence types and the number of matches in bestDict for each genome
            for sequenceType, matches in bestMatch[genome].iteritems():
                # If the number of matches for a profile matches the best number of matches
                if matches == sortedMatches:
                    # Iterate through the gene in the analysis
                    for gene in profileData[sequenceType]:
                        # Get the reference allele as above
                        refAllele = profileData[sequenceType][gene]
                        # As above get the reference allele split and ordered as necessary
                        if len(refAllele.split(" ")) > 1:
                            intRefAllele = map(int, refAllele.split(" "))
                            sortedRefAllele = " ".join(str(allele) for allele in sorted(intRefAllele))
                        else:
                            sortedRefAllele = refAllele
                        # Populate the sequence type dictionary with the genome, best match to profile, number of
                        # matches to the profile, gene, query allele(s), reference allele(s), and percent identity
                        MLSTseqType[genome][sequenceType][sortedMatches][gene][str(bestDict[genome][gene].keys()[0])][sortedRefAllele] = str(bestDict[genome][gene].values()[0])
            # Add the new profile to the profile file (if the option is enabled)
            if updateProfile:
                reProfilR(int(header), profileFile, geneList, genome)
        # Otherwise, the query profile matches the reference profile
        else:
            # Iterate through best match
            for sequenceType, matches in bestMatch[genome].iteritems():
                if matches == sortedMatches:
                    for gene in profileData[sequenceType]:
                        # Populate resultProfile with the genome, best match to profile, number of
                        # matches to the profile, gene, query allele(s), reference allele(s), and percent identity
                        resultProfile[genome][sequenceType][sortedMatches][gene][str(bestDict[genome][gene].keys()[0])] = str(bestDict[genome][gene].values()[0])
        dotter()


def reProfilR(numGenes, profileFile, geneList, genome):
    """Creates and appends new profiles as required"""
    global MLSTseqType
    global bestDict
    global resultProfile
    global profileData
    profileNumber = 0
    lastEntry = 0
    newProfile = ""
    # Iterate through MLSTseqType - it contains genomes with partial matches to current reference profiles
    # for genome in MLSTseqType:
    # Reset newProfile
    newProfile = ""
    # Find the last profile entry in the dictionary of profiles
    # Opens uses the command line tool 'tail' to look at the last line of the file (-1). This last line
    # is split on tabs, and only the first entry (the sequence type number) is captured
    profile = subprocess.check_output(['tail', '-1', profileFile]).split("\t")[0]
    # Split the _CFIA from the number - if there is no "_", the just use profile as the profile number
    try:
        profileNumber = int(profile.split("_")[0])
    except IndexError:
        profileNumber = int(profile)
    # If the number is less than 1000000, then new profiles have not previously been added
    if profileNumber < 1000000:
        # Set the new last entry number to be 1000000
        lastEntry = 1000000
    # If profiles have previously been added
    else:
        # Set last entry to the highest profile number plus one
        lastEntry = profileNumber + 1
    # As there can be multiple profiles in MLSTSeqType, this loop only needs to be performed once.
    seqCount = 0
    # Go through the sequence types
    for sequenceType in MLSTseqType[genome]:
        # Only do this once
        if seqCount == 0:
            # Set the newProfile string to start with the new profile name (e.g. 1000000_CFIA)
            newProfile = "%s_CFIA" % str(lastEntry)
            # The number of matches to the reference profile
            for numMatches in MLSTseqType[genome][sequenceType]:
                # The genes in geneList - should be in the correct order
                for gene in geneList:
                    # The allele for each gene in the query genome
                    for allele in MLSTseqType[genome][sequenceType][numMatches][gene]:
                        # Append the allele to newProfile
                        newProfile += "\t%s" % str(allele)
                        # Add the MLST results for the query genome as well as the new profile data to resultProfile
                        resultProfile[genome]["%s_CFIA" % str(lastEntry)][numGenes][gene][allele] = MLSTseqType[genome][sequenceType][numMatches][gene][allele].values()[0]
            seqCount += 1
    # Only perform the next loop if newProfile exists
    if newProfile:
        # Open the profile file to append
        appendFile = open(profileFile, "ab")
        # Append the newProfile to the end of the profile file
        appendFile.write("%s\n" % newProfile)
        # Close the profile file
        appendFile.close()
        # Remove the .json file with the old profile information
        jsonProfile = "%s.json" % os.path.splitext(profileFile)[0]
        try:
            os.remove(jsonProfile)
        except OSError:
            # If there are multiple new profiles, then this will be unable to delete the file each time; allow a pass
            pass
        # Re-run profilR with the updated files
        profilR(profileFile)


def blaster(path, cutoff, sequencePath, targetPath, organismPath, scheme, organism):
    """
    The blaster function is the stack manager of the module
    markers are the the target fasta folder that with be db'd and BLAST'd against strains folder
    out is the working directory where the blastxml folder will be placed
    name is the partial title of the csv output
    ALL PATHS REQUIRE TRAILING SLASHES!!!
    """
    # Time is used to calculate length of the analyses
    start = time.time()
    # Import global variables
    global count, genedict, blastpath, plusdict, updateProfile
    # Initialise genedict
    genedict = defaultdict(list)
    blastpath = []
    # Run organism chooser to allow the user to choose which databases to use
    # returns the organism name, lists of alleles and the profile
    alleles, profile, organismName, schemeName, profileFile = organismChooser(path, targetPath, organismPath, scheme, organism)
    print "[%s] Reading sequence profiles" % (time.strftime("%H:%M:%S"))
    profileDict, geneList = profilR(profileFile)
    # reset count to 0
    count = 0
    # Get the genome files into a list - note that they must be in the "sequences" subfolder of the path,
    # and the must have a file extension beginning with ".fa"
    strains = glob("%s*.fa*" % sequencePath)
    # # Create the threads for the BLAST analysis
    for i in range(len(strains)):
        threads = runblast(blastqueue)
        threads.setDaemon(True)
        threads.start()
    print "\n[%s] Creating necessary databases for BLAST" % (time.strftime("%H:%M:%S"))
    # Push targets to threads
    makedbthreads(alleles)
    print "\n[%s] Performing and parsing BLAST database searches" % (time.strftime("%H:%M:%S"))
    # reset count to 0
    count = 0
    # Make blastn threads and retrieve xml file locations
    blastnthreads(alleles, strains, cutoff)
    # reset count to 0
    count = 0
    print "\n[%s] Determining sequence types" % (time.strftime("%H:%M:%S"))
    # Determine sequence types
    sequenceTyper(profileDict, profileFile, geneList, updateProfile)
    # Parse the results into a report
    csvheader = ''
    # # Initialise variables
    row = ""
    rowcount = 0
    # Make the reports folder if necessary
    make_path("%sreports" % path)
    # Open the csv report - add the organism name, the scheme name and the date to keep reports unique
    with open("%sreports/%s_%s_results_%s.csv" % (path, organismName, schemeName, time.strftime("%Y.%m.%d.%H.%M.%S")), 'wb') as csvfile:
        # Get the header started
        csvheader = "Strain,SequenceType,Matches,"
        # Initialise the headerGenes variable
        headerGenes = ''
        for gene in geneList:
            # Append each gene to headerGenes
            headerGenes += "%s," % gene
        # Append headerGenes to the header
        csvheader += headerGenes
        # Write the header to the report
        csvfile.write(csvheader)
        # Iterate through all the query genomes
        for genome in plusdict:
            # Reset resultString to a newline
            resultString = "\n"
            # Reset the counts to 0
            sequenceTypeCount = 0
            resProfileCount = 0
            # If the genome is in resultProfile
            if genome in resultProfile:
                # If, as is the case for cgMLST, there are multiple identical profiles, then, the formatting of the
                # report is automatically changed
                if len(resultProfile[genome]) == 1:
                    # Iterate through the sequence types
                    for sequenceType in resultProfile[genome]:
                        # Put the genome and the sequence type in the result string
                        resultString += "%s,%s, " % (genome, sequenceType)
                        # Report the number of matches to the profile
                        for numMatches in resultProfile[genome][sequenceType]:
                            # Append these matches
                            resultString += "%s," % numMatches
                            # Go through the ordered list of genes
                            for gene in geneList:
                                # Add each allele to the result string
                                for allele in resultProfile[genome][sequenceType][numMatches][gene]:
                                    resultString += "%s," % allele
                    # Append the resultString to the report
                    csvfile.write(resultString)
                # If there are more than one identical profile
                else:
                    # Same as above, but only add the genome on the first line
                    for sequenceType in resultProfile[genome]:
                        # First sequence type
                        if resProfileCount == 0:
                            resultString += "%s,%s, " % (genome, sequenceType)
                            for numMatches in resultProfile[genome][sequenceType]:
                                resultString += "%s," % numMatches
                                for gene in geneList:
                                    for allele in resultProfile[genome][sequenceType][numMatches][gene]:
                                        resultString += "%s," % allele
                            resProfileCount += 1
                            csvfile.write(resultString)
                        # Subsequent sequence types
                        else:
                            resultString = "\n,"
                            resultString += "%s," % sequenceType
                            for numMatches in resultProfile[genome][sequenceType]:
                                resultString += "%s," % numMatches
                                for gene in geneList:
                                    for allele in resultProfile[genome][sequenceType][numMatches][gene]:
                                        resultString += "%s," % allele
                            resProfileCount += 1
                            csvfile.write(resultString)
            # The option to not update the profiles is optionally available, so do the same as above, but for genomes
            # that lack a perfect match to a reference profile
            # I'm not commenting here, as this is essentially identical to the code above. Ideally, I should find a
            # way to get this into a reusable function.
            else:
                if len(MLSTseqType[genome]) == 1:
                    for sequenceType in MLSTseqType[genome]:
                        resultString += "%s,%s," % (genome, sequenceType)
                        for numMatches in MLSTseqType[genome][sequenceType]:
                            resultString += "%s," % numMatches
                            for gene in geneList:
                                for allele in MLSTseqType[genome][sequenceType][numMatches][gene]:
                                    for refAllele in MLSTseqType[genome][sequenceType][numMatches][gene][allele]:
                                        if allele == refAllele:
                                            resultString += "%s," % allele
                                        # Since there are mismatches, show the expected allele in the reference profile
                                        else:
                                            resultString += "%s (%s)," % (allele, refAllele)
                    csvfile.write(resultString)
                else:
                    for sequenceType in MLSTseqType[genome]:
                        if sequenceTypeCount == 0:
                            resultString += "%s,%s," % (genome, sequenceType)
                            for numMatches in MLSTseqType[genome][sequenceType]:
                                resultString += "%s," % numMatches
                                for gene in geneList:
                                    for allele in MLSTseqType[genome][sequenceType][numMatches][gene]:
                                        for refAllele in MLSTseqType[genome][sequenceType][numMatches][gene][allele]:
                                            if allele == refAllele:
                                                resultString += "%s," % allele
                                            else:
                                                resultString += "%s (%s)," % (allele, refAllele)
                            sequenceTypeCount += 1
                            csvfile.write(resultString)
                        else:
                            resultString = "\n,"
                            resultString += "%s," % sequenceType
                            for numMatches in MLSTseqType[genome][sequenceType]:
                                resultString += "%s," % numMatches
                                for gene in geneList:
                                    for allele in MLSTseqType[genome][sequenceType][numMatches][gene]:
                                        for refAllele in MLSTseqType[genome][sequenceType][numMatches][gene][allele]:
                                            if allele == refAllele:
                                                resultString += "%s," % allele
                                            else:
                                                resultString += "%s (%s)," % (allele, refAllele)
                            sequenceTypeCount += 1
                            csvfile.write(resultString)

    # File cleanup
    tmpFiles = glob("%stmp/*" % path)
    for tmpFile in tmpFiles:
        try:
            os.remove(tmpFile)
        except OSError:
            raise
    # Calculate the elapsed time
    end = time.time() - start
    # Friendly exit statement
    print "\n[%s] Elapsed time for GeneSeeking is %.2f seconds with %.2f seconds per genome" \
          % (time.strftime("%H:%M:%S"), end, end/float(len(strains)))


# Run the blaster function
blaster(path, cutoff, sequencePath, allelePath, organismPath, scheme, organism)
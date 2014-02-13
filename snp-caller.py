import fileinput
from multiprocessing import Pool
import sys
import time
from math import log

fileprefix = 'any_prefix'

# constants
linesep     = '\n'
snpRate     = 0.001
ratioTiTv   = 2
bases       = {'A','C','G','T'}

numThreads = 4
linesToRead = 1000

pileupBases = bases.union({',','.'})


# builds dictionary of dictionary for priors
# {refbase: {read base: prior probability}}
def buildSimplePriorsDict(snpRate, ratioTiTv):
    bases = ['A','T','G','C']
    priors = {}
    for base1 in bases:
        priors[base1] = {}
        for base2 in bases:
            if base1 == base2:
                priors[base1][base2] = 1 - snpRate
            else:
                baseLst = [base1,base2]
                baseLst.sort()
                if baseLst == ['A','G'] or baseLst == ['C','T']: # transition
                    priors[base1][base2] = snpRate * ratioTiTv / (ratioTiTv + 1)
                else:
                    priors[base1][base2] = snpRate / (ratioTiTv + 1) / 2
    return priors

def calculateErrorProbs():
    errorProbs = {}
    for i in range(33,74):
        errorProbs[chr(i)] = 10 ** ((33-i) / 10.)
    return errorProbs


errorProbs  = calculateErrorProbs()
priors      = buildSimplePriorsDict(snpRate, ratioTiTv)


# this function calls the SNP for a single line of pileup and is passed to map
def process(line):

    data = line.split()

    # data[0]   [1] [2]     [3]         [4]         [5]
    # reference pos refBase coverage    read bases  baseQualities
    # chr22     1   G       1           ^~.         I
    refbase     = data[2].upper()
    qualities   = data[5]

    # initialize cumulative probabilities
    A = An = C = Cn = G = Gn = T = Tn = 1

    # iterate through bases in the reads
    readindex = 0
    for readbase in data[4].upper():

        # skip letters that are not bases, e.g. start or end of read markers
        if readbase in pileupBases:

            # skip bases with p(error)=1
            if qualities[readindex] == '!':
                readindex += 1
                continue

            # . or , means reference base
            if readbase == '.' or readbase == ',':
                readbase = refbase

            # look up error probability using quality score
            errorProb = errorProbs[qualities[readindex]]

            if readbase == 'A':
                A *= (1-errorProb)
                An *= errorProb
            elif readbase == 'C':
                C *= (1-errorProb)
                Cn *= errorProb
            elif readbase == 'G':
                G *= (1-errorProb)
                Gn *= errorProb
            else:
                T *= (1-errorProb)
                Tn *= errorProb

            readindex += 1

    probs = [
        ['A', priors[refbase]['A']*A*Cn*Gn*Tn],
        ['C', priors[refbase]['C']*An*C*Gn*Tn],
        ['G', priors[refbase]['G']*An*Cn*G*Tn],
        ['T', priors[refbase]['T']*An*Cn*Gn*T]]

    # consensus base is the one that maximizes the posterior probability
    [call,qual] = max(probs, key=lambda x: x[1])
    phred = int(-10 * log(qual,10))

    # print quasi-VCF format: CHROM  POS ID  REF ALT QUAL
    if call != refbase:
        return '\t'.join([data[0],data[1],'.',refbase,call,str(phred),'\n'])


# business end
with open(fileprefix+'.vcf','w') as output:
    output.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\n')

    queue = []
    pool = Pool(processes=numThreads)

    for line in fileinput.input():
        if len(queue) < linesToRead:
            queue.append(line)
        else:
            result = pool.map(process,queue)
            for item in result:
                if item:
                    output.write(item)
            queue = []

    for item in pool.map(process,queue):
        if item:
            output.write(item)
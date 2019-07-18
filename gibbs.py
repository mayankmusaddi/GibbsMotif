import random
import math
from weblogo import *

## a is a matrix containing starting index of motif in each sequence
## w is the size of the motif to be found

indices = { 'A' : 0, 'a' : 0, 'C' : 1, 'c' : 1 , 'T' : 2, 't' : 2, 'G' : 3, 'g' : 3}

## creating profile with sequence with index z ignored
def createProfile(sequences,a,w,z=-1):

	## initialiation of profile and background
	profile = [[0 for i in range(w)] for j in range(4)]
	background = [0]*4

	i = 0
	for sequence in sequences:
		if(i!=z):
			j = 0
			for base in sequence:
				if( j >= a[i] and j < a[i]+w):
					profile[ indices[base] ][ j-a[i] ] += 1
				else:
					background[indices[base]]+=1
				j+=1
		i+=1

	## normalize profile and background
	for i in range(0,len(profile)):
		for j in range(0,len(profile[0])):
			if(z!=-1):
				profile[i][j] = (profile[i][j]+0.5)/(len(sequences)-1+2)
			else:
				profile[i][j] = (profile[i][j]+0.5)/(len(sequences)+2)
	
	backgroundSum = sum(background)
	for i in range(0,len(background)):
		background[i] = (background[i]+0.5)/(backgroundSum+2)

	return profile,background

def sampleSeq(sequences,profile,background,w,z):
	Ax = []
	for i in range(0,len(sequences[z])-w):
		LR = 0
		for j in range(0,w):
			base = sequences[z][i+j]
			baseIn = indices[base.strip()]
			LR += math.log( (profile[baseIn][j]) / background[baseIn] )
		Ax.append(LR)
	return Ax

def getLikelihood(sequences,profile,background,a,w):
	L = 0
	i = 0
	for sequence in sequences:
		for j in range(0,len(sequence)):
			base = sequence[j]
			baseIn = indices[base.strip()]
			if( j>=a[i] and j<a[i]+w ):
				L += math.log(profile[baseIn][j-a[i]])
			else:
				L += math.log(background[baseIn])
		i+=1
	return L

# :
	## for reading fasta files
	def getSequences(fil):

	    lines = fil.readlines()
	    sequences = []

	    seqFlag = False
	    sequence = ''

	    for line in lines:
	        if line!="\n":
	            if line[0] == '>':
	                seqFlag = True
	                if sequence != '':
	                    sequences.append(sequence)
	                # print("SEQ ",sequence)
	                sequence = ''
	            else:
	                seqFlag = False
	            if seqFlag is False:
	                sequence += line[:-1]

	    if sequence != '':
	        sequences.append(sequence)
	    return sequences

	def localMinimaCheck(sequences,w,a,globalMaxAlignmentProb):

		buffer = 5
		na = a[:]
		for i in range(1,buffer):
			sa = a[:]
			for j in range(0,len(a)):
				if a[j]+i+w <= len(sequences[j]):
					sa[j] = a[j]+i

			profile,background = createProfile(sequences,sa,w,-1)
			AlignmentProb = getLikelihood(sequences,profile,background,sa,w)
			print("YO :: ",i,AlignmentProb)
			print("SA :: ",sa)

			if(AlignmentProb > globalMaxAlignmentProb):
				globalMaxAlignmentProb = AlignmentProb
				na = sa[:]
			else:
				break

		for i in range(1,buffer):
			sa = a[:]

			for j in range(0,len(a)):
				if a[j]-i >= 0:
					sa[j] = a[j]-i

			profile,background = createProfile(sequences,sa,w,-1)
			AlignmentProb = getLikelihood(sequences,profile,background,sa,w)
			print("YO :: ",-i,AlignmentProb)
			print("SA :: ",sa)


			if(AlignmentProb > globalMaxAlignmentProb):
				globalMaxAlignmentProb = AlignmentProb
				na = sa[:]
			else:
				break

		print("AFTER OPTIMIZATION")
		print(na)
		print(globalMaxAlignmentProb)

		i = 0
		for sequence in sequences:
			print(sequence[na[i]:na[i]+w])
			i+=1

def getMotif(sequences,w,MAXLOOP = 100,ITERATIONS = 500):

	globalMaxAlignmentProb = float('-inf')

	####Run Iterations with different random start position
	for i in range(0,ITERATIONS):

		####Initialize Random Alignment
		a = []
		for sequence in sequences:
			a.append(random.randint(0,len(sequence)-w-1))

		localMaxAlignmentProb = float('-inf')

		####Run Iterations with updates in start position
		for j in range(0,MAXLOOP):

			####Each sequence in the set has its motif start site updated
			for k in range(0,len(sequences)):
				profile,background = createProfile(sequences,a,w,k)
				Ax = sampleSeq(sequences,profile,background,w,k)
				a[k] = Ax.index(max(Ax))

			profile,background = createProfile(sequences,a,w,-1)

			AlignmentProb = getLikelihood(sequences,profile,background,a,w)

			if(AlignmentProb > localMaxAlignmentProb):
				localMaxAlignmentProb = AlignmentProb
			else:
				break

		if(localMaxAlignmentProb == globalMaxAlignmentProb):
			break
		elif(localMaxAlignmentProb > globalMaxAlignmentProb):
			globalMaxAlignmentProb = localMaxAlignmentProb

	i = 0
	motifs = []
	for sequence in sequences:
		motifs.append(sequence[a[i]:a[i]+w])
		i+=1

	return motifs,profile,globalMaxAlignmentProb

def plotMotif(motifs,title):
	seqMotifs = seq_io.array_io.read(motifs)
	seqMotifs.alphabet = Alphabet.which(seqMotifs)

	logodata = LogoData.from_seqs(seqMotifs)
	logooptions = LogoOptions(resolution = 300)
	logooptions.logo_title = title
	logoformat = LogoFormat(logodata, logooptions)

	motifJPG = logo_formatter.jpeg_formatter(logodata, logoformat)
	fil = open("./outputs/"+title+".jpg", "wb")
	fil.write(motifJPG)
	fil.close()

if __name__ == '__main__':
	FILENAME = 'MA0080.2.sites'
	
	fil=open("./sequences/"+FILENAME, "r")
	seqs = read_seq_data(fil)

	motifs,profile,acc = getMotif(seqs,7,1000,50000)

	print(acc)
	print(motifs)
	for pr in profile:
		print(pr)

	plotMotif(motifs,FILENAME)
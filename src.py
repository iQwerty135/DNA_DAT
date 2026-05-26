#ALL THE OUTPUTS ARE UPPERCASE

rna2amino = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y','UAC' : 'Y', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
"""dictionary rna2amino contains all the codon:amino-acid key:value pairs where:
        codon: a three-letter string object representing a codon
        amino-acid: an upper case one-letter string object representing the coded amino-acid
"""



def gcContent(dna):
    """Returns the GC content ratio of a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a real number between 0 and 1
    Example: gcContent("atcgttcaag") = 0.4
    """
    count = 0 
    if len(dna) == 0:
        return 0
    for base in dna.upper():
        if base in 'CG':
            count += 1
    return round(count / len(dna),1)   #ratio between the number of 'G','C' and the total number of bases 



def countCodon(dna, codon):
    """Returns the number of (non-overlapping) occurrences of a codon in a DNA sequence
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value: the integer number of instances of the target codon in dna
    Example: countCodon("aaaaaaaa", "aaa") = 2
    """
    dna , codon = dna.upper() , codon.upper()
    count = 0
    for i in range(0,len(dna),3):   #i: starting index of every codon
        if dna[i:i+3] == codon:
            count += 1
    return count



def countACG(dna):
    """Returns the number of nucleotides that are not T in a DNA sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: the integer number of nucleotides in dna that are not T
    Example: countACG("atcgttcaag") = 7
    """
    count = len(dna)
    for base in dna.upper():
        if base == 'T':
            count = count - 1
    return count



def printCodons(dna):
    """Prints the sequence of non-overlapping codons in dna, 
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printCodons("ggtacactgta") would print: ggt aca ctg
    """
    i = 0                       
    while i + 2 < len(dna):                     
        print(dna[i:i+3].upper(), end=' ')
        i += 3
    print()



def findCodon(dna, codon):
    """Returns the index of the first occurrence of codon in dna
    Parameters:
        dna: a string object representing a DNA sequence
        codon: a three-letter string object representing the codon to
                search for
    Return value:
        (if codon found): the integer index of the first occurrence of codon in dna
        (if codon not found): None
    Example: findCodon("ggtacactacgta", "tac") = 2 
    """
    index = dna.upper().find(codon.upper()) 
    if index != -1:
        return index
    return None


def findATG(dna):
    """Returns a list of all the positions of the codon ATG in dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a list of integer numbers
    Example: findATG("gatgtatgta") = [1,5] 
    """
    return findseq(dna,'ATG')



def printReadingFrames(dna):
    """Prints the sequences of non-overlapping codons in the dna
    with the three possible reading frames in separate columns
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: None
    Example: printReadingFrames("aggcctggc") should print
        agg    ggc    gcc
        cct    ctg    tgg
        ggc
    """
    if len(dna)>2:                                                                                   
        all_possible_codons = [dna[i:i+3].upper() for i in range(len(dna)-2)] + ['']*(2-len(dna)%3)  #add empty strings to avoid error raised by selection with index out of range
        for i in range(0,len(all_possible_codons),3):
            print(all_possible_codons[i],all_possible_codons[i+1],all_possible_codons[i+2],sep='\t')
    print()



def firstSSR(dna, seq):
    """Returns the length (number of repeats) of the first SSR in dna
    that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer number of seq repeats in the first SSR
        (if seq not found): 0
    Example: firstSSR("aggcctggcggcggc", "ggc") = 1
    """
    dna , seq = dna.upper() , seq.upper()
    index = dna.find(seq)                   #index of the first occurrence of seq in dna
    if index != -1 and seq != '':
        dna = dna[index+len(seq):]
        count = 1
        while seq == dna[:len(seq)]:        #iteration to count the number of repetition of the SSR
            dna = dna[len(seq):]
            count += 1
        return count
    return 0
    


def longestSSR(dna, seq):
    """Returns the length of the longest SSR in dna that repeats the sequence seq
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a short sequence of DNA
    Return value:
        (if seq found): the integer length of longest SSR in dna repeating seq
        (if seq not found): 0
    Example: longestSSR("aggcctggcggcggc", "ggc") = 3
    """
    index = findCodon(dna,seq)              #index of the first occurrence of seq in dna
    if index == None or seq == '':
        return 0
    len_firstSSR = firstSSR(dna, seq)
    return max(longestSSR(dna[index + (len_firstSSR * len(seq)):], seq), len_firstSSR)  #return the max value among all the SSR founded in a recursive way
    


def longestSSRdin(dna):
    """Finds the longest SSR in dna for all the possible dinucleotides
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: (if len(dna)>1): a pair (dins, len)
        dins: a list of two-letter string objects representing the dinucleotides with the longest SSR in dna
        len: integer representing the length of the longest SSR for all dincucleotides
                  (if len(dna)<2): ([], 0)
    Examples: longestSSRdin("ctctctgcgccacacaca") = (["ca"], 4)
              longestSSRdin("gagacacactctct") = (["ac", "ct"], 3)
    """
    if len(dna)<2:
        return ([],0)
    dna = dna.upper()
    all_din = [n1+n2 for n1 in 'ACTG' for n2 in 'ACTG']   #the dinucleotide at index 'i' in the list 'all_din'
    Longest = [longestSSR(dna,seq) for seq in all_din ]   #has the value of his longestSSR stored in 'Longest' at the same index                              
    maxlen = max(Longest)
    return ([all_din[i] for i in range(16) if Longest[i]==maxlen],maxlen) 

    

def complement(dna):
    """Returns the complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the complement of the DNA sequence
    Example: complement("acgtac") = "tgcatg"
    """
    complementary_dna = ''
    conjugated_bases = {'T':'A','A':'T','C':'G','G':'C'}
    for base in dna.upper():
        if base in conjugated_bases:
            complementary_dna += conjugated_bases[base]
        else :
            complementary_dna += 'N'   #it treats ambiguous symbols like clean
    return complementary_dna



def reverseComplement(dna):
    """Returns the reverse  complement of a dna sequence
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing the reverse complement of the DNA sequence
    Example: reverseComplement("acgtac") = "gtacgt"
    """
    return complement(dna)[::-1]    


    
def palindrome(dna):
    """Returns true if dna is the same as its reverse complement
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: Bool
    Example: palindrome("atat") = True
    """
    if len(dna)%2 == 0: #parity check 
        if clean(dna[0:len(dna)//2].upper()) == reverseComplement(dna[len(dna)//2:len(dna)]): #check half and half:clean and rev.comp. treat ambiguous symbol in the same way (subst with 'N')
            return True 
    return False



def dna2rna(dna):
    """Returns a copy of dna in which every "t" has been replaced by a "u"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: dna2rna("actgat") = "acugau"
    """
    rna = ''
    for base in dna.upper():
        if base != 'T':
            rna += base
        else:
            rna += 'U'
    return rna



def transcribe(dna):
    """Returns the RNA equivalent of the reverse complement of dna
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing an RNA sequence
    Example: transcribe("acgtac") = "guacgu"
    """
    return dna2rna(reverseComplement(dna))


    
def clean(dna):
    """Returns a new DNA string in which every character in dna
    that is not an "a", "c", "g", or "t" is replaced with an "n"
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object representing a clean DNA sequence
    Example: clean("goat") = "gnat"
    """
    c_dna = ''
    for base in dna.upper():
        if base in 'ATCG':
            c_dna += base
        else:
            c_dna += 'N'
    return c_dna



import random
Table={'A':('A'),'T':('T'),'C':('C'),'G':('G'),'*':('*'),'R':('A','G'),'Y':('C','T'),'K':('G','T'),'M':('A','C'),'S':('C','G'),'W':('A','T'),'B':('C','G','T'),'D':('A','G','T'),'H':('A','C','T'),'V':('A','C','G'),'N':('A','C','G','T')}
"""dictionary Table that contains all pairs key:value where:
    key : all the possible characters that can be found in dna
    value : all the possible nucleotide they represent
    special case :   '*':('*') to avoid error when is passed as a parameter a string that contains a character not present in Table as a key
"""
                        


def fix(dna):
    """Returns a DNA string in which each ambiguous symbol is replaced
    with one of the possible bases it represents, each with equal probability
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a string object representing a fixed DNA sequence
    Example: fix("arta") returns either "aata" or "agta"
        (with probability 1/2 for each one of these two possible returns)
    """
    if dna == '':
        return ''
    fixed = ''
    for base in dna.upper():
        if base not in Table:
            base = '*'
        fixed += random.choice(Table[base])            
    return fixed


def fixAll(dna):
    """Returns the list containing all possible DNA strings in which each
    ambiguous symbol in dna is replaced with all of the possible bases it
    represents; all combinations of replacements should be considered and
    returned within the list of fixed DNA strings
    Parameter:
        dna: a string object representing a DNA sequence with ambiguous simbols
    Return value: a list of string objects representing fixed DNA sequences
    Examples:   fixAll("arma") = ["aaaa", "aaca", "agaa", "agca"]
                fixAll("agata") = ["agata"]
    """
    if dna == '':
        return ['']
    ambiguous = dna[0].upper() if dna[0].upper() in Table else '*'
    return [substitute + remaining_cases for remaining_cases in fixAll(dna[1:]) for substitute in Table[ambiguous]]



def readFASTA(filename):
    """Reads a FASTA file and removes the header and all the newline characters
    Returns the DNA sequence contained in the file as a string
    Parameter:
        filename: the name of file containing a DNA sequence in the FASTA format
    Return value: a string object representing the DNA sequence in the file
    """
    dna = ''
    fasta_lines = open(filename).readlines()
    for line in fasta_lines[1:]:                #starts at index 1 to avoid the header  
        dna += line.strip()
    return dna



import urllib.request



def getFASTA(id):
    """Fetches the DNA sequence with the given id from the NCBI database
    and returns it as a string (header and newline characters must be removed)
    Parameter:
        id: a string object representing the identifier (NCBI accession number) of a DNA sequence
    Return value: a string object containing the dna sequence
    """
    url_file = urllib.request.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+id+'&rettype=fasta&retmode=text').readlines() 
    genome = ''
    for line in url_file[1:]:                   #starts at index 1 to avoid the header
        genome += line.decode('utf-8').strip()  
    return genome



def findseq(dna, seq):
    """Returns the list of indexes (starting positions) of all the
    occurrences of seq in dna
    Parameters:
        dna: a string object representing a DNA sequence
        seq: a string object representing a sequence of DNA
    Return value: a list of integer numbers
    """
    if seq == '':
        return []
    dna , seq = dna.upper() , seq.upper()
    L = [0]*len(seq)
    indexes = []
    j = 0 
    i = 1
    #algorithm to preprocess the pattern 
    while i < len(seq):
        if seq[j] == seq[i]:  
            j += 1
            L[i] = j
            i += 1
        elif j == 0:
            L[i] = 0
            i += 1
        else:
            j = L[j-1]
            
    i = 0
    j = 0
    #kmp search
    while i < len(dna):
        if dna[i] != seq[j]:
            if j == 0:
                i += 1
            else:
                j = L[j-1]
        else:
            i += 1
            j += 1
            if j == len(seq):
                indexes += [i-j]
                j = L[j-1]
    return indexes



def mark(dna):
    """Returns a new DNA string in which every start codon (atg) in dna
    is replaced with ">>>" and every stop codon (taa, tag, or tga) is
    replaced with "<<<"
    The function does not consider overlapping codons (but just the
    reading frame starting from offset 0)
    If two subsequent start codons are found (without a stop codon in
    between), only the first one should be replaced with ">>>", since
    codon atg also codes for the Methionine amino-acid
    Parameter:
        dna: a string object representing a DNA sequence
    Return value: a string object represnting a marked sequence of DNA
    Example: mark("ttgatggagatgcattagaag") = "ttg>>>gagatgcat<<<aag"
    """
    dna = dna.upper()
    marked_dna = ''
    marker = False          #for every iteration : True after finding a start codon and not yet a stop codon, False the opposite
    for e in range(0,len(dna),3):
        codon = dna[e:e+3]
        if not marker :
            if codon == 'ATG':
                codon = '>>>'
                marker = True
        elif marker :
            if codon in ['TGA','TAA','TAG']:
                codon = '<<<'
                marker = False
        marked_dna += codon
    return marked_dna



def proteins(marked_dna):
    """Returns the list of proteins traduced from marked_dna
    Proteins are represented as strings of amino-acids
    Proteins are obtained from the RNA sequences obtained from
    the dna enclosed between the markers ">>>" and "<<<"
    Parameter:
        marked_dna: a string object representing a marked sequence of DNA
    Return Value: a list of string objects representing proteins
    Example: proteins("ttg>>>gagcat<<<aagcag>>>aca>>>caccaacag<<<aga") = ["EH","HQQ"]
    """
    marked_dna = '>>><<<' + marked_dna.upper()                            #add '>>><<<' to avoid loss of genes after split('>>>')
    all_genes_with_stop = [gene for gene in marked_dna.split('>>>') if '<<<' in gene]   
    all_genes = [dna2rna(gene.split('<<<')[0]) for gene in all_genes_with_stop]     #elimine the sequence after the first '<<<'
    all_genes = [[rna2amino[gene[i:i+3]] for i in range(0,len(gene),3) if gene[i:i+3] in rna2amino] for gene in all_genes]  
    all_genes = [''.join(gene) for gene in all_genes if not gene==[]]
    return all_genes


""" The following lines define the base settings and functions for plotting through the turtle module
You do not need to add any code to the functions plot and bar!
"""

width = 1200		# width of the window
cols = width // 6	# number of columns of text
height = 600		# height of the window
rows = height // 100	# number of rows of text

import turtle

def plot(tortoise, index, value, window):
	"""Plots GC fraction value for window ending at position index."""
	
	if (index == window) or (index - window + 1) // cols != (index - window) // cols:
		tortoise.up()	
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)
		tortoise.down()
	else:
		tortoise.goto((index - window + 1) % cols, \
		              (index - window + 1) // cols + 0.7 + value * 0.25)

		
def bar(tortoise, index, rf):
	"""Draws a colored bar over codon starting at position index in
	   reading frame rf. Puts the turtle's pen up and down to
	   handle line breaks properly."""
	   
	tortoise.up()
	tortoise.goto(index % cols, index // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 1) % cols, (index + 1) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)
	tortoise.up()
	tortoise.goto((index + 2) % cols, (index + 2) // cols + (rf + 1) / 5)
	tortoise.down()
	tortoise.forward(1)


"""
Complete the following two functions to accomplish the required tasks
"""

def orf(dna, rf, tortoise):
    """Finds and draws all ORFs in the reading frame rf
    Blue bars begin and end on start and stop codons
    Parameters:
        dna: a string object representing a sequence of DNA
        rf: reading frame offset (it's value can be either 0, 1, or 2)
        tortoise: the drawing turtle
    Return value: None
    """
    marker = False    #True: when we find a start codon, False: after placing the blue bar over a stop codon
    dna = mark(dna[rf:])
    tortoise.pencolor('red')
    for e in range(0,len(dna)-2,3):
        if not marker :     
            if dna[e] == '>':
                tortoise.pencolor('blue')
                marker = True
            bar(tortoise,e+rf,rf)
        elif marker :
            bar(tortoise,e+rf,rf)
            if dna[e] == '<':
                marker = False
                tortoise.pencolor('red')    	



def gcFreq(dna, window, tortoise):
    """Computes and plots the GC frequency in dna over a sliding window
    Parameters:
        dna: a string object representing a sequence of DNA
        window: integer size of the sliding window
        tortoise: the drawing turtle
    Return value: None
    """

    # draws red lines at 0.5 above the sequence:
	
    tortoise.pencolor('red')
    for index in range(len(dna) // cols + 1):
    	tortoise.up()
    	tortoise.goto(0, index + 0.825)
    	tortoise.down()
    	if index < len(dna) // cols:
    		tortoise.goto(cols - 1, index + 0.825)
    	else:
    		tortoise.goto((len(dna) - window) % cols, index + 0.825)
    tortoise.up()
    tortoise.pencolor('blue')


    #our code below
    
    for e in range(len(dna)-window+1):  #stop plotting the gc content when the number of nucleotides in dna[e:e+window] is smaller than the window
        plot(tortoise,e+window-1,gcContent(dna[e:e+window]),window)



"""
The viewer functuion calls the functions you have coded to build the final plot!
"""

def viewer(dna):
    """Displays GC content and ORFs in 3 forward reading frames."""

    dna = dna.upper()   # makes everything upper case

    tortoise = turtle.Turtle()
    screen = tortoise.getscreen()
    screen.setup(width, height)     # makes a long, thin window
    screen.setworldcoordinates(0, 0, cols, rows) # scales coord system so 1 char fits at each point
    screen.tracer(100)
    tortoise.hideturtle()
    tortoise.speed(0)
    tortoise.up()

    # Prints the DNA string in the window:
	
    for index in range(len(dna)):
    	tortoise.goto(index % cols, index // cols)
    	tortoise.write(dna[index], font = ('Courier', 9, 'normal'))
		
    # Finds ORFs in forward reading frames 0, 1, 2:
	
    tortoise.width(5)
    for rf in range(3):
    	orf(dna, rf, tortoise)
		
    # Plots GC frequency:
	
    tortoise.width(1)
    gcFreq(dna, 5, tortoise)

    screen.update()

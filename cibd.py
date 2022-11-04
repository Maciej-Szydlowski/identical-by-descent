#!/usr/bin/python3
import subprocess
import os
import itertools




#This algorithm finds the most common chromosomal segments
#shared identical-by-descent (IBD). It was used to find common 
#segments for the paper by Szydlowski: 'A clue to the etiology of
#disorders of sex development from identity-by-descent analysis
#in dogs with cryptic relatedness' published in Animal Genetics.

#To search for most common chromosomal segment (of IBD2 type:
#both alleles shared identical-by-descent), the algorithm
#considers all subsets of individuals, starting from the largest
#subsets of size n-1, n-2, n-3 etc. For example, in the full set
#of n=11 individuals, the algorithm searches for common segments 
#shared by all 11 individuals and then generates 11 subsets of 10
#individuals, 55 subsets of 9 individuals, 165 subsets of 8 individuals,
#and so on. To find chromosomal segments shared identical-by-descent
#in a set of individuals, the shared fragments (of IBD2 type) are first
#detected in all pairwise comparisons using KING (version 2.2.7).

#You need to have BEDTOOLS installed on your computer.
#https://bedtools.readthedocs.io/en/latest/

#Before using the algorithm, please run KING (version 2.2.7). 
#You need the output file 'king.segments'.
#Reference: Manichaikul et al. 2010, Bioinformatics 26: 2867-2873.
#Now specify the output from KING
king_segments_file = 'example.king.segments'

#Specify IDs of individuals to consider.  
#In this example I consider six individuals and all subsets of them:
SAMPLES = [ '6260', '6311', '6419', '6466', '6798', '6826' ]

#Specify minimum subset size.
minimum_subset_size = 3

#The output file 'cibd.out.segments.bed' gives:
#Chrom, Start, End, Segment size, IDS of individuals sharing the segment (joined with '&').

#The output file 'cibd.out.sets.txt' gives:
#Segment size, Subset size, IDs in the subset.





def findsubsets(s, n):
  return list(itertools.combinations(s, n))


def doit( cmd ):
    p = subprocess.Popen( cmd, shell=True )
    os.waitpid( p.pid, 0 )


def getSEG( infile, outfile, sample1, sample2, ibdtype ):    
    fout = open( outfile, 'w+' )
    with open( infile ) as fin :
        next( fin )
        i = 0 
        for line in fin :
            line = line.strip() 
            (FID1, ID1, FID2, ID2, IBDType, Chr, StartMB, StopMB, StartSNP, StopSNP, N_SNP, Length ) = line.split( '\t' )
            if IBDType == ibdtype :
                start = int( float(StartMB) * 1000000 )
                stop  = int( float(StopMB ) * 1000000 )
                if( (ID1==sample1 and ID2==sample2) or (ID1==sample2 and ID2==sample1) ) :
                    d = stop - start + 1    
                    tag = sample1+'_'+sample2+'_'+IBDType
                    string = Chr + '\t' + str( start ) + '\t' + str( stop ) + '\t' + str( d ) + '\t' + tag + '\n'
                    fout.write( string )
    fout.close()
 
def size( infile ):
    size = 0
    with open( infile ) as fin :
        for line in fin :
            line = line.strip() 
            vv = line.split(  )
            chrom=vv[0]; start=vv[1]; end=vv[2];
            len = int( end ) - int( start ) + 1
            size = size + len
    return size ;         





fout1  = open( 'cibd.out.sets.txt', 'w+' ) 
fout2  = open( 'cibd.out.segments.bed', 'w' )
for k in range( len( SAMPLES ), minimum_subset_size - 1, -1 ) :     
     sets = findsubsets( SAMPLES, k )
     for aset in sets :
         nsample = len( aset )
         for i in range( nsample-1 ):
             j = i + 1
             getSEG( king_segments_file, 'temp1', aset[i], aset[j], 'IBD2' )
             if i == 0:
                 cmd = 'cp temp1 product.bed' ; doit( cmd )
             if i >  0 :
                 cmd = 'sort -k1,1 -k2,2n temp1 > temp1s' ; doit( cmd )
                 cmd = 'sort -k1,1 -k2,2n product.bed > temp2s' ; doit( cmd )
                 cmd = 'bedtools intersect -a temp1s -b temp2s > temp3' ; doit( cmd )
                 cmd = 'cp temp3 product.bed' ; doit( cmd )
         w = size( 'product.bed' )   
         if w > -1 :
             st = ''
             for ii in range( len(SAMPLES) ) :
                 znak = '----'
                 if SAMPLES[ii] in aset :
                     znak = SAMPLES[ii]
                 st = st + znak + '\t'
             st = str(w) + '\t' + str(nsample) + '\t' + st
             print(st) 
             if w > 0 :
               tag = '&'.join( aset );             
               with open( 'product.bed' ) as fin :
                   for line in fin :
                       line = line.strip() 
                       ( chrom, start, end, nic, tagnic ) = line.split('\t')
                       dl = int(end) - int(start); dl=str(dl)
                       wers = '\t'.join( [chrom, start, end, dl, tag] )
                       #print( wers )
                       wers = wers + '\n'
                       fout2.write( wers )
               st = st + '\n'
               fout1.write( st )

fout1.close()
fout2.close()

cmd = 'sort -k1,1nr  cibd.out.sets.txt > cibd.out.sets.sorted.txt' ; doit( cmd )
cmd = 'sort -k1,1 -k2,2n cibd.out.segments.bed > cibd.out.segments.sorted.bed'; doit( cmd )
cmd = 'rm temp1 temp1s temp2s temp3 product.bed'; doit( cmd )

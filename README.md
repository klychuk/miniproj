
# miniproject A
#E. coli is a model organism. 
#The strain E. coli K-12 was isolated from the stool of a convalescent diphtheria patient in 1922, and it has been used in the lab for nearly 100 years. 
#It was one of the first organisms to have its whole genome sequenced.
#It was such a big deal it was published in 1997 in Science! 
#Most researchers either buy their strain from a stock collection or know the history of their strain (
#These strains are going to Inevitably evolve over time. 
#Recently, researchers have returned back to these K-12 strains to resequence them and 
#ones that were derived (evolved) from the original K-12 strain


# This code is comparing Ecoli strains to the original K-12 strain
#It uses wget, prokka, spades, bowtie, tophat and cufflinks

# Method accepts:
#SRR1 is the strain in question in the format 'SRR0000000'
#RefSRR is the K12 strain in the same format 
#fqn is the # of fastq files
#user is the username of the user i.e. 'klychuk'
#genus is the escherichia

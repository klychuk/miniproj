import os
import csv

path = os.getcwd()
os.mdkir(path+'/optionA_Karson_Lychuk')
os.chdir()

#SRR1 is the strain in question in the format 'SRR0000000'
#RefSRR is the K12 strain in the same format 
#fqn is the # of fastq files
#user is the username of the user i.e. 'klychuk'
#genus is the escherichia
def mini_proj(SRR1,Ref_SRR,fqn,user,genus):
   #various formatting needed to construct the wget command
    SRRR = SRR1 + ".sra"
    SRR = SRR1 + '.fasta'
    g = 'wget'
    seq_first_three = SRR1[3] + SRR1[4] + SRR1[5]
   
    #supplies the wget command for the strain in question SRR1 aka 8185310
    basicform = " ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/SRR/SRR"
    seq_get = g + basicform + seq_first_three + '/' + SRR1 + "/" + SRRR
   
    #uses the os to run the command made
    os.system(seq_get)
   
    #fastq dump command to split and or convert the sra file
    #converts the sra of the strain to a fastq file
    seq_fastq_d = "fastq-dump -I --split-files " + SRRR
    os.system(seq_fastq_d)
    with open(OptionA.log, 'w') as output:
        output.write(seq_get + '\n' + seq_fastq_d)
        
    #sorts by the number of fastq files present to specify the process
    #in our case it is one so it should be -s
    fqn = 1
    if fqn == 1:
        di = "-s " + SRR
    elif fqn == 2:
        di = "-1 " + SRR + '_1 ' + '-2 ' + SRR + '_2'
   
    
    # the spades command, using the algorithm, proper output dir and input which is our fasta file
    # WHAT SPADES DOES LMAO
    spades = 'spades -k 55,77,99,127 -t 2 --only-assembler ' + di + " -o home/" + user + '/spades_' + SRR1

    #print(spades)
    os.system(spades)
    
    #to interprety our spades findings a bit of python code is needed
    #first the contigs file is opened and parsed
    #then the # of contigs with a len > 1000 is counted and appended to a list
    #then the total basepairs in the list of contigs is counted for a final bp count
    seqs = []
    real = []
    tot = 0
    directory = "/home/" + user 
    contig_info = directory + '/contigs.fasta'
    
    with open(contig_info, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seqs.append(str(record.seq))
            for s in seqs:
                if len(s) > 1000:
                    real.append(s)
                    #print(len(real))    
                    for s in real:
                        tot += len(s)
                        #print(tot)
    a1 = '\n There are ' + len(real) +' contigs > 1000 in the assembly.\n'
    a2 = 'There are 883431035'+ tot+' bp in the assembly.'
    with open(OptionA.log, 'a') as output:
        output.write(spades + '\n' + a1 + a1)
       
    
    #assembles the prokka command with the outfir, specificed genus and location of the contig file
    #WHAT IT DOES LMAO
    proka = 'prokka --outdir prok --prefix ' + genus + ' home/' + user + '/contigs.fasta\n'
    os.system(proka)
    
    #view the stats of the prokka file and compare to the refseq data
  with open('/home/klychuk/mp/mp/escherichia.txt') as f:
    found = 0
    dif = 0
    content = f.readlines()
    #print(content)
    #print(content[0])
    for i in range(len(content)):
        c = content[i]
        #print(c[0])
        if c[0] == 'C':
            icds = c[5:9]
            cds = int(icds)
            ognum = 4140
            found = cds - ognum
        if c[1] == 'R':
            trna = c[6:8]
            itrna = int(trna)
            #print(itrna)
            tfind = 89
            dif = tfind - itrna

    a3  = 'Prokka found found an additional %d CDS and %d less tRNA than RefSeq.\n'%(found,dif)
    with open(OptionA.log, 'a') as output:
        output.write(proka + '\n' + a3)
    
    SRR2 = 'SRR1411276'
    ref_first_3 = SRR2[3] + SRR2[4] + SRR2[5] 
    refget = g + basicform + ref_first_3 + '/' + SRR2 + '/' + SRR2 + '/' + SRR2 + '.sra\n'
    os.system(refget)
    
    os.system(refget)
    ref_fastqd = "fastq-dump -I -spit-files " + SRR2 + '.sra'
    os.system(ref_fastqd)
    
    with open(OptionA.log, 'a') as output:
        output.write(refget + '\n' + ref_fastqd)
    
    #in order to run tophat bowtie to map
    #creates a bowtie index for our ref seq
    fnaf = 'NC_000913.fna'
    b2bcom = 'bowtie2-build ' + fnaf + ' EcoliK12' 
    os.system(b2bcom)
    
    #run bowtie with the RNA seq
    #assembles the reference genome
    bt= 'bowtie2 -x Ecolik12 -U '
    btcom = bt + SRR2 + '.fastq -S EcoliK12.sam'
    os.system(btcom)
    
    #now using the bowtie output run cufflinks
    #analyzeds the expression quantification between the strain and the ref
    clc = 'cufflinks -p 2 -accepted_hits.bam'
    os.system(clc)
    
     with open(OptionA.log, 'a') as output:
        output.write(brbcom + '\n' + btcom + '\n' + clc)
        
    #gathers the fpkm from isoforms
    with open('isoforms.fpkm_tracking', mode = 'r') as _file:
        csv_reader = csv.DictReader(_file, delimiter='\t')
        lc = 0
        for row in csv_reader:
            if lc == 0:
                lc += 1
           # print(row.keys())
           else:
               a4 = (f'{row["locus"]}\t {row["FPKM"]}')

    with open(OptionA.log, 'a') as output:
        output.write(a4)
                                  
     
    
#mini_proj('SRR8185310','SRR1411276',1,'klychuk','escherichia')

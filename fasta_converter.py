import re
filename='/data/siting/lncRNAseq-longsee-20180208/DGElncRNA_Stie/siD2vsNC/sub_test.fa'
savename='/data/siting/lncRNAseq-longsee-20180208/DGElncRNA_Stie/siD2vsNC/reformat_sub.fa'

fasta=open(filename, 'r' )
save=open(savename, 'w')

for line in fasta:
	if line: # to make sure you aren't at an empty line
		save.write(re.sub( r'^(>[^\s]+).*', r'\1', line)
                save.close()
                fasta.close()

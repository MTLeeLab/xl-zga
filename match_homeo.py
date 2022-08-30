def match(infile):
	f=open(infile, 'r')
	for line in f:
		fields=line.strip().split('\t')
		lgenesU=fields[4].split(',')
		sgenesU=fields[9].split(',')
		lgenesD=fields[5].split(',')
		sgenesD=fields[10].split(',')
		countUS=0
		countUS2=0
		countDS=0
		countDS2=0
		for gene in lgenesU:
			if gene in sgenesU:
				countUS=countUS+1
			else:
				if gene in sgenesD:
					countUS2=countUS2+1
		for gene in lgenesD:
			if gene in sgenesD:
				countDS=countDS+1
			else:
				if gene in sgenesU:
					countDS2=countDS2+1
		lcoord=fields[1]+':'+fields[2]+'-'+fields[3]
		scoord=fields[6]+':'+fields[7]+'-'+fields[8]
		print(fields[0],lcoord,scoord,countUS+countUS2,countDS+countDS2,sep='\t')

match('landmark/agg/both_on_total.txt')

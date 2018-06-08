def gen_haplotypes(builder, genotype, haplotypes):
	# Base Case: weve reached end of the string, so add to list
	if len(genotype) == 0:
		haplotypes.append(builder);
	# inductive case:
	else:
		if genotype[0] == '0':
			# Always 0, add to builder and recurse
			builder += '0'
			gen_haplotypes(builder, genotype[1:], haplotypes)
		elif genotype[0] == '1': 
			# Recurse with adding both a 0 and a 1
			gen_haplotypes(builder+'1', genotype[1:], haplotypes)
			gen_haplotypes(builder+'0', genotype[1:], haplotypes)
		elif genotype[0] == '2': 
			# Always 1, add to builder and recurse
			builder += '1'
			gen_haplotypes(builder, genotype[1:], haplotypes)


# Returns all phases for a specific genome
def gen_phases(genotype):
	haplotypes = []
	gen_haplotypes('', genotype, haplotypes)
	phases = []
	# Special case: only 1 haplotype, so create a phase of that haplotype with itself
	if len(haplotypes) == 1:
		phz = (haplotypes[0], haplotypes[0])
		phases.append(phz)
	else: 
		# Match 0 with n, 1 with n-1, etc
		n = int(len(haplotypes)/2)
		for i in range(0, n):
			h1 = haplotypes[i]
			h2 = haplotypes[-(i+1)]

			phases.append((h1, h2))
	return phases
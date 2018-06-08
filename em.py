import pandas 
import sys
from tqdm import tqdm
from generate_phases import gen_phases

def main():
	#1: parse args
	if len(sys.argv) == 2:
		data_set = sys.argv[1]
	else:
		print("Invalid args. Usage: py em.py [path_to_data_set]")
		return -1

	#2: load data into a pandas dataframe
	#	rows: SNPs
	#	cols: individual

	# numerical column names for pandas dataframe
	_CSV_COLUMNS = [
		'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
		'20', '21', '22', '23', '24', '25', '26', '27', '28', '29',
		'30', '31', '32', '33', '34', '35', '36', '37', '38', '39',
		'40', '41', '42', '43', '44', '45', '46', '47', '48', '49'
	]
	df = pandas.read_csv(data_set, delimiter=' ', names=_CSV_COLUMNS )

	#3: Generate subset of genotypes
	SUBSET_SIZE = 16
	NUM_SNPS = df.shape[0]
	NUM_INDIVIDUALS = df.shape[1]
	NUM_SUBSETS = int(NUM_SNPS / SUBSET_SIZE) 

	# For each subset:
	for i in tqdm(range(0, NUM_SUBSETS+1)):
		SNP_BEGIN = i*SUBSET_SIZE
		SNP_END = SNP_BEGIN + SUBSET_SIZE

		#3.1: Generate genotypes for each individual
		# genotypes: [a list of lists] for all genotypes in subset range
		genotypes = []
		# For each individual:
		for j in range(0, NUM_INDIVIDUALS):
			genotype = df.iloc[SNP_BEGIN:SNP_END, j].tolist()
			genotypes.append(''.join(str(x) for x in genotype))


		#3.2: Generate a genotype and haplotype dictionary 
		# create 2 dictionary data structures:
		#	genotype# -> ([(pha,ses)], probabilities)
		#	haplotype -> [p1, p2]
		genotype_dict = {}  # note, if mapping genotypes directly, there will be duplicates
							# unaccounted for. instead, just map their index
		haplotype_dict = {}

		# For each genotype:
		for j in range(0, len(genotypes)):
			#from generate_phases.py
			phases = gen_phases(genotypes[j])
			# create probabilities list to correspond with phases
			probabilities = []
			num_phases = len(phases)
			for k in range(0, num_phases):
				probabilities.append(1.0/num_phases)

			# add to genotype dictionary
			genotype_dict[j] = (phases, probabilities)

		# iterate through the genotype dictionary and update the haplotype dict
		# for each genotype
		for geno in genotype_dict:
			#for each phase in phase_list
			phase_pos = 0
			for phz in genotype_dict[geno][0]:
				#first haplotype in pair:
				hap0 = phz[0]
				if hap0 not in haplotype_dict:
					#init haplotype in dict
					haplotype_dict[hap0] = [0.0, 0.0] 

				#second haplotype in pair:
				hap1 = phz[1]
				if hap1 not in haplotype_dict:
					#init haplotype in dict
					haplotype_dict[hap1] = [0.0, 0.0] 

				#increment phase pos
				phase_pos += 1


		# all haplotypes are found, update probabilities 
		init_haplotype_probability = 1.0 / len(haplotype_dict)
		for haplo in haplotype_dict:
			haplotype_dict[haplo][0] = init_haplotype_probability


		#3.3: Run iterative EM, updating probabilities for phases and haplotypes on each iteration
		num_iters = 10

		# For each iter:
		for j in range(0, num_iters):
			# For each genotype:
			for geno in genotype_dict:
				total_p = 0
				# For each phase:
				phase_pos = 0
				for phz in genotype_dict[geno][0]:
					hap0 = phz[0]
					hap1 = phz[1]
					p = haplotype_dict[hap0][0] * haplotype_dict[hap1][0]
					total_p += p
					genotype_dict[geno][1][phase_pos] = p #update probability
					#increment phase pos
					phase_pos += 1

				# For each phase:
				phase_pos = 0
				for phz in genotype_dict[geno][0]:
					hap0 = phz[0]
					hap1 = phz[1]
					p = genotype_dict[geno][1][phase_pos] / total_p
					genotype_dict[geno][1][phase_pos] = p
					haplotype_dict[hap0][1] += p
					haplotype_dict[hap1][1] += p
					#increment phase pos
					phase_pos += 1

			for haplo in haplotype_dict:
				new_p = haplotype_dict[haplo][1] / (2*len(genotype_dict))
				haplotype_dict[haplo][0] = new_p
				haplotype_dict[haplo][1] = 0

		#3.4: Find phase with highest probability and print it to file line by line
		# For each genotype:
		pm_list = []
		for geno in genotype_dict:
			p_max = 0.0
			phase_max = ('', '')
			# For each phase:
			phase_pos = 0
			for phz in genotype_dict[geno][0]:
				if genotype_dict[geno][1][phase_pos] >= p_max:
					phase_max = genotype_dict[geno][0][phase_pos]
					p_max = genotype_dict[geno][1][phase_pos]
				#increment phase pos
				phase_pos += 1

			# max phase for genotype found
			pm_list.append(phase_max)

		# We found all max phase, now write it to file
		#50 genotypes, 100 lines
		num_lines = SUBSET_SIZE
		if i == NUM_SUBSETS:
			#we are on the last set, so there are going to be less
			snps_left = NUM_SNPS - i*SUBSET_SIZE
			num_lines = snps_left

		with open('output.txt', 'a') as the_file:
			for j in range(0, num_lines):
				line = ''
				for phz in pm_list:
					line += phz[0][j] + ' ' + phz[1][j] + ' '
				line += '\n'
				the_file.write(line)


	#FIN
	print ('Finished processing file: %s...' % data_set)



if __name__ == "__main__":
	main()
import pysam
import pandas as pd

# Define your list of individuals
list_ind = [
    # ... your list of individual ids]
]

def calculate_heterozygosity_from_vcf(vcf_file, list_ind):
    vcf = pysam.VariantFile(vcf_file)
    # Check if individuals in list are in VCF
    for ind in list_ind:
        if str(ind) not in vcf.header.samples:
            raise ValueError(f"Individual {ind} not found in VCF")

    heterozygosity = []
    for ind in list_ind:
        hets = 0
        total_sites = 0
        for record in vcf.fetch():
            if record.samples[ind]["GT"]:
                genotype = record.samples[ind]["GT"]
                if genotype == (0, 1) or genotype == (1, 0):  # Heterozygous
                    hets += 1
            total_sites += 1
        heterozygosity.append(hets / total_sites)

    # Create a DataFrame with results
    df = pd.DataFrame({
        "Individual": list_ind,
        "heterozygosity": heterozygosity
    })
    return df

def calculate_allele_frequencies(vcf):
    frequencies = {}
    for record in vcf.fetch():
        ref_count_A = sum(1 for sample in parental_population_A if record.samples[sample]["GT"] in [(0,0)])
        alt_count_A = sum(1 for sample in parental_population_A if record.samples[sample]["GT"] in [(1,1)])
        total_A = ref_count_A + alt_count_A
        
        ref_count_B = sum(1 for sample in parental_population_B if record.samples[sample]["GT"] in [(0,0)])
        alt_count_B = sum(1 for sample in parental_population_B if record.samples[sample]["GT"] in [(1,1)])
        total_B = ref_count_B + alt_count_B
        
        if total_A > 0 and total_B > 0:
            freq_A = ref_count_A / total_A
            freq_B = ref_count_B / total_B
         
            frequencies[record.pos] = (freq_A, freq_B)
    
    return frequencies

def calculate_hybrid_index_for_individual(individual, frequencies, vcf):
    values = []
    for record in vcf.fetch():
        if record.pos not in frequencies:
            continue
        freq_A, freq_B = frequencies[record.pos]
        genotype = record.samples[individual]["GT"]
        
        if genotype == (0, 0):
            value = freq_A / (freq_A + freq_B)
        elif genotype == (1, 1):
            value = freq_B / (freq_A + freq_B)
        else:
            value = 0.5  # heterozygous
            
        values.append(value)

    return sum(values) / len(values)

# Load VCF
vcf = pysam.VariantFile("try_output.vcf.gz")

parental_population_A = ["9990", "9992", "9982", "9984", "10146", "9956"]
parental_population_B = ["5075", "8606", "8770"]
hybrid_population = list_ind  # assuming the hybrid population is the same as list_ind

frequencies = calculate_allele_frequencies(vcf)

hybrid_indices = {}
for individual in hybrid_population:
    hybrid_index = calculate_hybrid_index_for_individual(individual, frequencies, vcf)
    hybrid_indices[individual] = hybrid_index

# Now, calculate heterozygosity
vcf_file_path = "output_fixed_all.vcf.gz"
results_df = calculate_heterozygosity_from_vcf(vcf_file_path, list_ind)
results_df['hybrid_index'] = results_df['Individual'].map(hybrid_indices)

# Save the combined results to a CSV file
results_df.to_csv("combined_results.csv", sep=",", index=False)

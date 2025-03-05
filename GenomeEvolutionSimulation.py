#!/usr/bin/env python
# coding: utf-8

# In[5]:


pip install msprime simuPOP biopython numpy pandas matplotlib seaborn


# In[1]:


import numpy as np

# Population parameters
POPULATION_SIZE = 100
GENOME_LENGTH = 1000

# Create a homogeneous genome ('A' at all positions)
population = ["A" * GENOME_LENGTH for _ in range(POPULATION_SIZE)]

# Display the first 100 bases of the first individual
print(population[0][:100])


# In[2]:


import random

MUTATION_RATE = 0.001  # Probability of mutation for each genome position
NUCLEOTIDES = ["A", "T", "C", "G"]

def introduce_mutations(genome):
    genome_list = list(genome)
    for i in range(len(genome_list)):
        if random.random() < MUTATION_RATE:
            genome_list[i] = random.choice([n for n in NUCLEOTIDES if n != genome_list[i]])
    return "".join(genome_list)

# Mutate the entire population
population = [introduce_mutations(genome) for genome in population]

# Display the first 100 bases after mutations
print(population[0][:100])


# In[3]:


def recombine(parent1, parent2):
    crossover_point = random.randint(0, GENOME_LENGTH)
    child = parent1[:crossover_point] + parent2[crossover_point:]
    return child

# Create a new population through crossover
new_population = []
for _ in range(POPULATION_SIZE):
    parent1, parent2 = random.sample(population, 2)
    child = recombine(parent1, parent2)
    new_population.append(child)

population = new_population


# In[4]:


NUM_GENERATIONS = 100  # Liczba pokoleń

for generation in range(NUM_GENERATIONS):
    population = [introduce_mutations(genome) for genome in population]  # Mutacje
    new_population = []
    
    for _ in range(POPULATION_SIZE):
        parent1, parent2 = random.sample(population, 2)
        child = recombine(parent1, parent2)
        new_population.append(child)
    
    population = new_population

    if generation % 10 == 0:
        print(f"Generacja {generation} - przykładowy genom: {population[0][:50]}")


# In[5]:


from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# Function to calculate the frequency of each nucleotide in the population
def compute_nucleotide_frequencies(population):
    positions = list(zip(*population))  # Transpose to analyze positions in the genome
    frequencies = [{nucleotide: positions[i].count(nucleotide) for nucleotide in NUCLEOTIDES} for i in range(GENOME_LENGTH)]
    return frequencies

# Get nucleotide frequencies
frequencies = compute_nucleotide_frequencies(population)

# Visualization of genetic diversity
plt.figure(figsize=(12, 6))
sns.heatmap([[freq[nuc] for nuc in NUCLEOTIDES] for freq in frequencies], cmap="viridis", xticklabels=NUCLEOTIDES)
plt.xlabel("Nucleotide")
plt.ylabel("Genome Position")
plt.title("Genetic Diversity After Simulation")
plt.show()


# In[ ]:


from Bio import Entrez, SeqIO

Entrez.email = "zuzanna.lubas@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="NC_000852", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

print(f"Actual genome: {record.seq[:100]}")


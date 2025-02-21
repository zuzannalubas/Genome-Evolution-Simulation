#!/usr/bin/env python
# coding: utf-8

# In[5]:


pip install msprime simuPOP biopython numpy pandas matplotlib seaborn


# In[6]:


import numpy as np

# Parametry populacji
POPULATION_SIZE = 100
GENOME_LENGTH = 1000

# Tworzymy jednorodny genom ('A' na wszystkich pozycjach)
population = ["A" * GENOME_LENGTH for _ in range(POPULATION_SIZE)]

# Wyświetlamy pierwsze 100 zasad pierwszego osobnika
print(population[0][:100])


# In[7]:


import random

MUTATION_RATE = 0.001  # Prawdopodobieństwo mutacji na każdą pozycję genomu
NUCLEOTIDES = ["A", "T", "C", "G"]

def introduce_mutations(genome):
    genome_list = list(genome)
    for i in range(len(genome_list)):
        if random.random() < MUTATION_RATE:
            genome_list[i] = random.choice([n for n in NUCLEOTIDES if n != genome_list[i]])
    return "".join(genome_list)

# Mutujemy całą populację
population = [introduce_mutations(genome) for genome in population]

# Wyświetlamy pierwsze 100 zasad po mutacjach
print(population[0][:100])


# In[8]:


def recombine(parent1, parent2):
    crossover_point = random.randint(0, GENOME_LENGTH)
    child = parent1[:crossover_point] + parent2[crossover_point:]
    return child

# Tworzymy nową populację przez krzyżowanie
new_population = []
for _ in range(POPULATION_SIZE):
    parent1, parent2 = random.sample(population, 2)
    child = recombine(parent1, parent2)
    new_population.append(child)

population = new_population


# In[9]:


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


# In[12]:


from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

# Funkcja obliczająca częstość występowania każdej zasady w populacji
def compute_nucleotide_frequencies(population):
    positions = list(zip(*population))  # Transpozycja, aby analizować pozycje w genomie
    frequencies = [{nucleotide: positions[i].count(nucleotide) for nucleotide in NUCLEOTIDES} for i in range(GENOME_LENGTH)]
    return frequencies

# Pobieramy częstości nukleotydów
frequencies = compute_nucleotide_frequencies(population)

# Wizualizacja różnorodności genetycznej
plt.figure(figsize=(12, 6))
sns.heatmap([[freq[nuc] for nuc in NUCLEOTIDES] for freq in frequencies], cmap="viridis", xticklabels=NUCLEOTIDES)
plt.xlabel("Nukleotyd")
plt.ylabel("Pozycja w genomie")
plt.title("Różnorodność genetyczna po symulacji")
plt.show()


# In[13]:


from Bio import Entrez, SeqIO

Entrez.email = "zuzanna.lubas@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="NC_000852", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
handle.close()

print(f"Rzeczywisty genom: {record.seq[:100]}")


# In[ ]:





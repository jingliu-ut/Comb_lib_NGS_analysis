# Variants Analysis

Here, we determine the number of unique dna and protein sequences in each
library given a range of abundance filters.


```python
import os
import glob
import pandas as pd
import dask.dataframe as dd
```

Get root directory


```python
root_dir = "/media/scratch/post_analysis"
print(f"Root directory: {root_dir}")
```

    Root directory: /media/scratch/post_analysis


Set `home_folder`


```python
home_folder = "merged/QCend"
home_dir = os.path.join(root_dir, home_folder)
print(f"Home directory: {home_dir}")

```

    Home directory: /media/scratch/post_analysis/merged/QCend


Set directory based on `folder` and `output_folder` specified below.


```python
folder = "variant_analysis"
output_folder = "variant_analysis"

os.chdir('/media/scratch/post_analysis/')
os.chdir(os.path.join(home_dir, folder))
work_dir = os.getcwd()
print(f"Current working directory: {work_dir}")
output_dir = os.path.join(home_dir, output_folder)
print(f"Output directory: {output_dir}")

```

    Current working directory: /media/scratch/post_analysis/merged/QCend/variant_analysis
    Output directory: /media/scratch/post_analysis/merged/QCend/variant_analysis


Get files by specifying `name`.


```python
name = "*.parquet"

def get_files(name):
    print("Getting files")
    files = sorted([file for file in glob.glob(name)])
    for file in files:
        print(file)
    return files


files = get_files(name)
print("Done")

```

    Getting files
    non_conserved_Ga.parquet
    non_conserved_H3.parquet
    Done


Read `parquet` files into DataFrame.


```python
for file in files:
    if "Ga" in file:
        print(f"Getting {file} for Ga")
        Ga = file
    elif "H3" in file:
        print(f"Getting {file} for H3")
        H3 = file

print(f"Reading Ga_dna")
Ga_dna = dd.read_parquet(Ga, engine='pyarrow')['seq'].compute()

print(f"Reading H3_dna")
H3_dna = dd.read_parquet(H3, engine='pyarrow')['seq'].compute()

print("Done")

```

    Getting non_conserved_Ga.parquet for Ga
    Getting non_conserved_H3.parquet for H3
    Reading Ga_dna
    Reading H3_dna
    Done



```python
Ga_dna

```




    0         TTTCCTTGTGTT
    1         CATACGGCGTGG
    2         TATGGGCCTATG
    3         ATTGCTGCGCGG
    4         GGGGGGGTGTAG
                  ...     
    554322    AGTAATAATGCG
    554323    TTGTGTTATGAT
    554324    AGTGAGTGGAGT
    554325    CATATTTCTATT
    554326    CGTGCTTGTGTT
    Name: seq, Length: 47705940, dtype: object




```python
H3_dna

```




    0         TAGTTTAGTTTT
    1         TGTTGTTTTCTT
    2         TTTCCGGATTAG
    3         GATAATTTTAAT
    4         TGTCGTAGGTAT
                  ...     
    542951    ATTGAGGAGTGT
    542952    GAGTATAGTAAG
    542953    GGGGGTATGGGG
    542954    TCGTTTCTGGGT
    542955    TTTAATATGTCT
    Name: seq, Length: 65779138, dtype: object



Translate dna sequence to protein sequence.


```python
def translate(seq):
    aa = ""
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        aa += table[codon]
    return aa

# print("Converting Ga_dna to Ga_aa")
# Ga_aa = Ga_dna.apply(translate)
# print("Converting H3_dna to H3_aa")
# H3_aa = H3_dna.apply(translate)
# print("Done")
```


```python
# Ga_aa
```


```python
# H3_aa
```

Applying `nnk_filter` to filter for sequences that conform to NNK
randomization at the designated positions for mutagenesis.


```python
def nnk_filter(df):
    seqs = df.to_list()
    filter = [seq for seq in seqs if not ("A" in seq[2::3] or "C" in seq[2::3])]
    filtered_df = pd.DataFrame()
    filtered_df['seq'] = filter
    return filtered_df

print("Applying nnk_filter to Ga_dna")
nnk_Ga_dna = nnk_filter(Ga_dna)['seq']
print("Applying nnk_filter to H3_dna")
nnk_H3_dna = nnk_filter(H3_dna)['seq']
print("Done")
```

    Applying nnk_filter to Ga_dna
    Applying nnk_filter to H3_dna
    Done



```python
# print("Convert nnk_Ga_dna to nnk_Ga_aa")
# nnk_Ga_aa = nnk_Ga_dna.apply(translate)
# print("Convert nnk_H3_dna to nnk_H3_aa")
# nnk_H3_aa = nnk_H3_dna.apply(translate)
```

    Done


Count number of unique dna sequences for unfiltered dna sequences.


```python
print(f"Counting unique dna variants for Ga_dna")
Ga_count_dna = Ga_dna.value_counts()
print(f"Counting unique dna variants for H3_dna")
H3_count_dna = H3_dna.value_counts()

# print(f"Counting unique aa variants for Ga_aa")
# Ga_count_aa = Ga_aa.value_counts()
# print(f"Counting unique aa variants for H3_aa")
# H3_count_aa = H3_aa.value_counts()
print("Done")

```

    Counting unique dna variants for Ga_dna
    Counting unique dna variants for H3_dna
    Done


Count number of unique dna sequences for nnk-filtered dna sequences.


```python
print(f"Counting unique dna variants for nnk_Ga_dna")
nnk_Ga_count_dna = nnk_Ga_dna.value_counts()
print(f"Counting unique dna variants for nnk_H3_dna")
nnk_H3_count_dna = nnk_H3_dna.value_counts()

# print(f"Counting unique aa variants for nnk_Ga_aa")
# nnk_Ga_count_aa = nnk_Ga_aa.value_counts()
# print(f"Counting unique aa variants for nnk_H3_aa")
# nnk_H3_count_aa = nnk_H3_aa.value_counts()
print("Done")
```

    Counting unique dna variants for nnk_Ga_dna
    Counting unique dna variants for nnk_H3_dna
    Done


For a range of abundance filters, the dna sequence calls are translated into
protein sequences. The number of unique dna and protein sequences are then
reported for each library.

NNK filtered versus non-NNK filtered numbers and proportions (%) of dna and
protein sequences call are also reported.



```python
def count_variants(num):
    for i in range(num):
        abundance = i + 1
        Ga_dna_data = Ga_count_dna[Ga_count_dna >= abundance]
        H3_dna_data = H3_count_dna[H3_count_dna >= abundance]
        Ga_aa_data = convert2aa(Ga_dna_data, "Ga")
        H3_aa_data = convert2aa(H3_dna_data, "H3")
        nnk_Ga_dna_data = nnk_Ga_count_dna[nnk_Ga_count_dna >= abundance]
        nnk_H3_dna_data = nnk_H3_count_dna[nnk_H3_count_dna >= abundance]
        nnk_Ga_aa_data = convert2aa(nnk_Ga_dna_data, "Ga")
        nnk_H3_aa_data = convert2aa(nnk_H3_dna_data, "H3")
        data = [Ga_dna_data, Ga_aa_data,
                H3_dna_data, H3_aa_data,
                nnk_Ga_dna_data, nnk_Ga_aa_data,
                nnk_H3_dna_data, nnk_H3_aa_data]
        names = ["Ga_dna_data", "Ga_aa_data",
                "H3_dna_data", "H3_aa_data",
                "nnk_Ga_dna_data", "nnk_Ga_aa_data",
                "nnk_H3_dna_data", "nnk_H3_aa_data"]
        report(data, names, abundance)
    return None

def convert2aa(s, lib):
    dna_s = s.index.to_series()
    filter = dna_s.values
    if lib == "Ga":
        dna_df = Ga_dna.to_frame()
    elif lib == "H3":
        dna_df = H3_dna.to_frame()
    filtered_dna_df = dna_df[dna_df['seq'].isin(filter)]
    aa_s = filtered_dna_df['seq'].apply(translate)
    aa_count = aa_s.value_counts()
    return aa_count

def report(data, names, abundance):
    print(f"\nAbundance filter = {abundance}")
    for i in range(len(data)):
        d = data[i]
        name = names[i]
        unique = len(d)
        print(f"Number of unique variants in {name} = {unique}")
        coverage = None
        if "dna" in name:
            coverage = (unique / (32**4)) * 100
        elif "aa" in name:
            coverage = (unique / (21**4)) * 100
        print(f"Coverage for {name} = {coverage}")
    return None

count_variants(5)
print("\nDone")

```

    
    Abundance filter = 1
    Number of unique variants in Ga_dna_data = 1169738
    Coverage for Ga_dna_data = 111.5549087524414
    Number of unique variants in Ga_aa_data = 194481
    Coverage for Ga_aa_data = 100.0
    Number of unique variants in H3_dna_data = 1144589
    Coverage for H3_dna_data = 109.15651321411133
    Number of unique variants in H3_aa_data = 194418
    Coverage for H3_aa_data = 99.96760609005507
    Number of unique variants in nnk_Ga_dna_data = 1044226
    Coverage for nnk_Ga_dna_data = 99.58515167236328
    Number of unique variants in nnk_Ga_aa_data = 194481
    Coverage for nnk_Ga_aa_data = 100.0
    Number of unique variants in nnk_H3_dna_data = 1047797
    Coverage for nnk_H3_dna_data = 99.92570877075195
    Number of unique variants in nnk_H3_aa_data = 194417
    Coverage for nnk_H3_aa_data = 99.96709190100832
    
    Abundance filter = 2
    Number of unique variants in Ga_dna_data = 1059245
    Coverage for Ga_dna_data = 101.01747512817383
    Number of unique variants in Ga_aa_data = 194481
    Coverage for Ga_aa_data = 100.0
    Number of unique variants in H3_dna_data = 1062319
    Coverage for H3_dna_data = 101.31063461303711
    Number of unique variants in H3_aa_data = 194308
    Coverage for H3_aa_data = 99.91104529491312
    Number of unique variants in nnk_Ga_dna_data = 1041737
    Coverage for nnk_Ga_dna_data = 99.34778213500977
    Number of unique variants in nnk_Ga_aa_data = 194481
    Coverage for nnk_Ga_aa_data = 100.0
    Number of unique variants in nnk_H3_dna_data = 1046038
    Coverage for nnk_H3_dna_data = 99.7579574584961
    Number of unique variants in nnk_H3_aa_data = 194308
    Coverage for nnk_H3_aa_data = 99.91104529491312
    
    Abundance filter = 3
    Number of unique variants in Ga_dna_data = 1043171
    Coverage for Ga_dna_data = 99.48453903198242
    Number of unique variants in Ga_aa_data = 194481
    Coverage for Ga_aa_data = 100.0
    Number of unique variants in H3_dna_data = 1046560
    Coverage for H3_dna_data = 99.8077392578125
    Number of unique variants in H3_aa_data = 194145
    Coverage for H3_aa_data = 99.8272324802937
    Number of unique variants in nnk_Ga_dna_data = 1039821
    Coverage for nnk_Ga_dna_data = 99.16505813598633
    Number of unique variants in nnk_Ga_aa_data = 194481
    Coverage for nnk_Ga_aa_data = 100.0
    Number of unique variants in nnk_H3_dna_data = 1043106
    Coverage for nnk_H3_dna_data = 99.47834014892578
    Number of unique variants in nnk_H3_aa_data = 194145
    Coverage for nnk_H3_aa_data = 99.8272324802937
    
    Abundance filter = 4
    Number of unique variants in Ga_dna_data = 1038653
    Coverage for Ga_dna_data = 99.05366897583008
    Number of unique variants in Ga_aa_data = 194481
    Coverage for Ga_aa_data = 100.0
    Number of unique variants in H3_dna_data = 1039648
    Coverage for H3_dna_data = 99.1485595703125
    Number of unique variants in H3_aa_data = 193907
    Coverage for H3_aa_data = 99.7048554871684
    Number of unique variants in nnk_Ga_dna_data = 1038002
    Coverage for nnk_Ga_dna_data = 98.99158477783203
    Number of unique variants in nnk_Ga_aa_data = 194481
    Coverage for nnk_Ga_aa_data = 100.0
    Number of unique variants in nnk_H3_dna_data = 1038884
    Coverage for nnk_H3_dna_data = 99.07569885253906
    Number of unique variants in nnk_H3_aa_data = 193907
    Coverage for nnk_H3_aa_data = 99.7048554871684
    
    Abundance filter = 5
    Number of unique variants in Ga_dna_data = 1036015
    Coverage for Ga_dna_data = 98.80208969116211
    Number of unique variants in Ga_aa_data = 194478
    Coverage for Ga_aa_data = 99.99845743285977
    Number of unique variants in H3_dna_data = 1033243
    Coverage for H3_dna_data = 98.5377311706543
    Number of unique variants in H3_aa_data = 193503
    Coverage for H3_aa_data = 99.49712311228346
    Number of unique variants in nnk_Ga_dna_data = 1035864
    Coverage for nnk_Ga_dna_data = 98.78768920898438
    Number of unique variants in nnk_Ga_aa_data = 194478
    Coverage for nnk_Ga_aa_data = 99.99845743285977
    Number of unique variants in nnk_H3_dna_data = 1033055
    Coverage for nnk_H3_dna_data = 98.51980209350586
    Number of unique variants in nnk_H3_aa_data = 193503
    Coverage for nnk_H3_aa_data = 99.49712311228346
    
    Done



```python
import os
print('h')
print(os.getcwd())

```

    h
    /media/scratch/post_analysis/post_analysis/post_analysis


# Data conversion for further processing

FASTA files corresponding to each variant library (Ga and H3) are converted
into csv files to enable easy separation of bases at non-conserved (mutated)
and conserved base positions. After separation, the mutated and non-mutated
bases are written into parquet file format (uses RAPIDS and/or DASK on NVIDIA
GPUs) to enable much faster read and write speeds compared to csv format.
The mutated bases are then concatenated into a single file.




```python
import os
import glob
import csv
import dask.dataframe as dd
```

Setting `folder` and `output_folder`.


```python
root_dir = "/media/scratch/post_analysis"
print(f"Root directory: {root_dir}")
```

Set `home_folder`


```python
home_folder = "merged/QCend"
home_dir = os.path.join(root_dir, home_folder)
print(f"Home directory: {home_dir}")

```

Set directory based on `folder` and `output_folder` specified below.


```python
folder = "data/untrimmed"
output_folder = "data/untrimmed"

os.chdir('/media/scratch/post_analysis/')
os.chdir(os.path.join(home_dir, folder))
work_dir = os.getcwd()
print(f"Current working directory: {work_dir}")
output_dir = os.path.join(home_dir, output_folder)
print(f"Output directory: {output_dir}")
```

Get files by specifying `name`.


```python
name = "*.fasta"

def get_files(name):
    print("Getting files")
    files = sorted([file for file in glob.glob(name)])
    for file in files:
        print(file)
    return files


files = get_files(name)
print("Done")

```

Load fasta files and convert them to csv files.



```python
def convert_fasta(file):
    print(f"Retrieving information from {file}")
    f = open(file, 'r')
    out = file.split('.')
    out[1] = "csv"
    out = '.'.join(out)
    lines = f.read().splitlines()
    f.close()

    seqs = [line for line in lines if not line.startswith('>')]

    print("Splitting sequences into nucleotides")
    nucleotides = [list(seq) for seq in seqs]

    # Write output
    os.chdir(output_dir)
    print(f"Writing sequence dataframe to {out}")
    with open(out, 'w') as f:
        w = csv.writer(f)
        w.writerows(nucleotides)
    return None


for file in files:
    os.chdir(work_dir)
    # Convert fasta file to csv file
    convert_fasta(file)
print("Done")

```

Set directory based on `folder` and `output_folder` specified below.


```python
folder = "data/untrimmed"
output_folder = "data/trimmed"

os.chdir('/media/scratch/post_analysis/')
os.chdir(os.path.join(home_dir, folder))
work_dir = os.getcwd()
print(f"Current working directory: {work_dir}")
output_dir = os.path.join(home_dir, output_folder)
print(f"Output directory: {output_dir}")
```

Get files by specifying `name`.


```python
name = "*.csv"

def get_files(name):
    print("Getting files")
    files = sorted([file for file in glob.glob(name)])
    for file in files:
        print(file)
    return files


files = get_files(name)
print("Done")
```

Remove non-conserved (mutated) region and export sequence data to parquet files.
Need to manually specify `index` in `trim_df`.


```python
# QCend files
def trim_df(file):
    print(f"Loading dataframe from {file}")
    df = dd.read_csv(file, header=None)

    if "Ga" in file:
        index = [36, 37, 38, 39, 40, 41, 42, 43, 44, 48, 49, 50]
    elif "H3" in file:
        index = [10, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, 27]

    print("Trimming non-conserved region")
    trimmed_df = df.drop(df.columns[index], axis=1)
    return trimmed_df
```


```python
def write_parquet(df):
    df.columns = [str(i) for i in range(df.columns.size)]
    out = file.split('.')
    out[1] = "parquet"
    out = '.'.join(out)
    out = "trimmed_" + out
    print(f"Writing dataframe to {out}")

    os.chdir(output_dir)
    df.to_parquet(out,
                  engine='pyarrow',
                  write_metadata_file=False)


for file in files:
    os.chdir(work_dir)
    df = trim_df(file)
    write_parquet(df)
print("Done")

```

Set directory based on `folder` and `output_folder` specified below.



```python
folder = "data/untrimmed"
output_folder = "variant_analysis"

os.chdir('/media/scratch/post_analysis/')
os.chdir(os.path.join(home_dir, folder))
work_dir = os.getcwd()
print(f"Current working directory: {work_dir}")
output_dir = os.path.join(home_dir, output_folder)
print(f"Output directory: {output_dir}")
```

Get files by specifying `name`


```python
name = "Ga*.csv"

def get_files(name):
    print("Getting files")
    files = sorted([file for file in glob.glob(name)])
    for file in files:
        print(file)
    return files


files = get_files(name)
print("Done")
```

Writing non-conserved (mutated) region into parquet files.


```python
def non_conserved(file):
    print(f"Loading dataframe from {file}")
    trim_df = dd.read_csv(file, header=None)

    if "Ga" in file:
        index = [36, 37, 38, 39, 40, 41, 42, 43, 44, 48, 49, 50]
    elif "H3" in file:
        index = [10, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, 27]

    df = trim_df.iloc[:, index]
    return df

dfs = []
for file in files:
    df = non_conserved(file)
    df.columns = [str(i) for i in range(df.columns.size)]
    dfs.append(df)

result = dd.concat(dfs)


```


```python

variants = result.apply(''.join, axis=1, meta=('seq', 'string'))

```


```python
def write_parquet(df):
    out = "non_conserved_H3.parquet"
    print(f"Writing dataframe to {out}")

    os.chdir(output_dir)
    df.to_parquet(out,
                  engine='pyarrow',
                  write_metadata_file=False)

write_parquet(variants.to_frame())

```

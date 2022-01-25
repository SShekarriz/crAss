import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], delimiter = "\t")

df = df[df['contig_type'].str.contains('Complete')]

contigs= df["contig"]
print(contigs)


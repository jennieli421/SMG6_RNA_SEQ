import os
import sys
import pandas as pd

def add_sample_type(row):
	if row['Genotype'] == "Smg6 Flox/Flox; RosaCreER tg/+" and row['TREATMENT'] == "no treatment":
		return 'control'
	elif row['Genotype'] == "Smg6 Flox/Flox; RosaCreER +/+" and row['TREATMENT'] == "4-OHT":
		return 'control-4-OHT'
	elif row['Genotype'] == "Smg6 Flox/Flox; RosaCreER tg/+" and row['TREATMENT'] == "4-OHT":
		return 'Smg6-iKO'

"""
Add sample type column if the metadata doesn't have 
"""
def sample_type(metadata, output):

	md = pd.read_csv(metadata)

	md['Samplename'] = md.apply(lambda row: add_sample_type(row), axis=1)

	#only keep columns 'Run' and 'Samplename'
	md = md[['Run', 'Samplename']]

	md.to_csv(output, sep=' ', index=False, header=False) 


if __name__ == '__main__':
	# two arguments are: ${metadata} ${output}
	sample_type(sys.argv[1], sys.argv[2]) 

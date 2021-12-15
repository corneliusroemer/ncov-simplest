#%%
import pandas as pd
import numpy as np
from fuzzywuzzy import process
#%%
country_counts = pd.read_csv("tmp/country_counts.tsv",sep="\t")
country_counts.sort_values(by='count',ascending=False)
# %%
country_pop = pd.read_csv("tmp/pop.csv",header=None,names= ["country","pop"])
#%%
def map_country(country):
    return process.extractOne(country,country_pop['country'])[0]

# %%
country_counts['join_key'] = country_counts.country.apply(map_country)

# %%
df = country_counts.join(country_pop.set_index('country'),on='join_key')
df.set_index('country',inplace=True)
# %%
df['seqpercapita'] = df['count'] / df['pop']
df
# %%
df.sort_values(by='seqpercapita',ascending=False)

# %%
# Rule: at most 10 seq per million
# Calculate max
# Calculate min of actual and max
# Sum

df['max'] = np.ceil(40 * df['pop'] / 10e6)
df['actual'] = df[['max','count']].min(axis=1)
df['overflow'] = df['actual'] - df['count']
df.overflow
# %%
df.actual.sum()

# %%
df['max'].astype('int').to_csv("tmp/max.tsv",sep="\t")

max_df = df
# %%
# Load metadata
# Group by country
# Randomly select up to max from each group
df = pd.read_csv("builds/metadata.tsv",sep="\t")
#%%
df_weights = pd.read_csv("tmp/max.tsv",sep="\t").set_index('country')
df_weights['sequenced'] = df.country.value_counts()
df_weights['chosen'] = df_weights[['max','sequenced']].min(axis=1)
df_weights
# %%
out = df.groupby('country').apply(
    lambda x: x.sample(
        n=df_weights.loc[x.name,'chosen'],
        replace=False
        )
    )

# %%

#%%
from typing_extensions import Required
import pandas as pd
import numpy as np
from fuzzywuzzy import process
import click
#%%

@click.command()
@click.option('-i', '--input-metadata', required=True, type=click.Path(exists=True), help='Path to the input metadata file.')
@click.option('-o','--output-metadata', required=True, type=click.Path(exists=False), help='Path to the output metadata file.')
@click.option('-p', '--population', required=True, type=click.Path(exists=True), help='Path to the population file.')
def subsample(input_metadata, output_metadata, population):
    """Subsample a metadata file proportional to population"""
    # %%
    # input_metadata = "builds/21K/metadata.tsv"
    meta = pd.read_csv(input_metadata,sep="\t",low_memory=False)
    meta['clock_deviation'] = pd.to_numeric(meta['clock_deviation'], errors='coerce')
    meta.dropna(subset=['clock_deviation'], inplace=True)
    meta = meta[meta.clock_deviation < 10]
    meta = meta[meta['QC_rare_mutations'] == 'good']
    meta = meta[meta['QC_snp_clusters'] == 'good']
    meta = meta[meta['QC_missing_data'] == 'good']
    # %%
    seq_per_country = meta.groupby('country').virus.count()
    seq_per_country.sort_values(ascending=False,inplace=True)
    seq_per_country = seq_per_country.reset_index()
    seq_per_country['rewritten'] = seq_per_country.country
    seq_per_country.loc[seq_per_country.country == 'USA', 'rewritten'] = 'United States'
    seq_per_country.rename(columns={'virus':'count'},inplace=True)
    seq_per_country
    # %%
    country_pop = pd.read_csv(population,header=None,names= ["country","pop"])
    country_pop
    #%%
    def map_country(country):
        return process.extractOne(country,country_pop['country'])[0]
    # %%
    seq_per_country['join_key'] = seq_per_country.rewritten.apply(map_country)
    seq_per_country
    # %%
    df = seq_per_country.join(country_pop.set_index('country'),on='join_key')
    df.set_index('country',inplace=True)
    df
    # %%
    df['seqpercapita'] = df['count'] / df['pop']
    df
    # %%
    df.sort_values(by='seqpercapita',ascending=False)

    # %%
    ratio = 10e-5
    while True:
        df['max'] = np.ceil(ratio * df['pop'])
        df['actual'] = df[['max','count']].min(axis=1).astype(int)
        if df['actual'].sum() < 1500:
            print(f"{df['actual'].sum()} sequences")
            print(f"{ratio} sequences per population")
            break
        ratio *= 0.9
    df

    # %%
    out = meta.groupby('country').apply(
        lambda x: x.sample(
            n=df.loc[x.name,'actual'],
            replace=False
            )
        ).to_csv(output_metadata,sep="\t",index=False)


if __name__ == '__main__':
    subsample()
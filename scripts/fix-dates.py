#%%
import pandas as pd
from treetime.utils import numeric_date
from datetime import datetime as dt
from datetime import timedelta
#%%
def danish_daterange(datestring):
    """
    Converts a datestring in the format YYYY-MM-DD
    to numeric date range plus 7 days
    """
    start_date = dt.strptime(datestring, '%Y-%m-%d')
    end_date = start_date + timedelta(days=7)
    return f"[{numeric_date(start_date)}:{numeric_date(end_date)}]"

def add_xx_to_incomplete_date(datestring):
    """
    Adds xx to incomplete dates
    So: 2023-08 becomes 2023-08-XX
    """
    if len(datestring) == 7:
        return datestring + "-XX"
    else:
        return datestring
#%%
# Read in all the metadata as strings in pandas
df = pd.read_csv("builds/BA.2.86/metadata.tsv",sep="\t",low_memory=False,dtype=str)
#%%
# Create treetime date column from date column
df["raw_date"] = df.date
#%%
# If country is "Denmark" apply danish_daterange to date
df.loc[df.country == "Denmark", "date"] = df.loc[df.country == "Denmark", "date"].apply(danish_daterange)
# %%
# Add xx to incomplete dates
df.date = df.date.apply(add_xx_to_incomplete_date)

df.to_csv("builds/BA.2.86/treetime_metadata.tsv",sep="\t",index=False)

# %%

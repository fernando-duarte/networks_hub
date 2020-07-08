import wrds
import sys
import re
import pandas as pd

# Change to match your WRDS username and WRDS cloud file structure
wrds_username = ''
cwd = ''

db = wrds.Connection()

# sql_query = 'select cusip, permco from CRSPM.STOCKNAMES'
# query_result = db.raw_sql(sql_query)
query_result = db.get_table(library='crspa', \
                            table = 'stocknames', columns = ['cusip', 'permco'])

query_result.to_csv(cwd + 'permco_cusip.csv')

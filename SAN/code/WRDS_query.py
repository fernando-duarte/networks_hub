import wrds
import sys
import re
import pandas as pd


wrds_username = ''
cwd = '/home/frb-ny/' + wrds_username + '/'

db = wrds.Connection()

sql_query = 'select cusip, permco from CRSPM.STOCKNAMES'
query_result = db.raw_sql(sql_query)

query_result.to_csv(cwd + 'permco_cusip.csv')
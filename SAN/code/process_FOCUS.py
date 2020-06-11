##############################################
#
# 	FOCUS Raw File Processing
#	Collin Jones
#	Purpose: Take in SIFMA data on aggregated balance sheet info for brokers and dealers. Combine disparate data
#		files, produce series of interest for Glasserman upper bound.
#
################################################


import pandas as pd 
import re
import datetime
import numpy as np


all_df = pd.DataFrame()
for year in [2012, 2013, 2014, 2015, 2016, 2017]:
	for quarter in [1, 2, 3, 4]:
		for firms in ['Top 10', 'Top 11-25']:

			#year = 2012
			#quarter = 1
			df = pd.read_excel('../input/SIFMA ' + str(year)[-2:] + ' Q' + str(quarter) + '.xls', sheetname=firms, skiprows=[0])[['BALANCE SHEET', 'Unnamed: 4']]
			df.loc[[('TOTAL LIABILITIES' in str(d)) and ('&' not in str(d)) for d in df['BALANCE SHEET']], 'side'] = 'Liabilities'
			df.loc[['TOTAL ASSETS' in str(d) for d in df['BALANCE SHEET']], 'side'] = 'Assets'
			df['side'] = df['side'].bfill()
			df = df.dropna(subset = ['Unnamed: 4', 'BALANCE SHEET', 'side'])

			maincat = df['BALANCE SHEET'].apply(lambda x: not bool(re.match(r' *-', x)))
			df.loc[maincat, 'main'] = df['BALANCE SHEET']
			df['main'] = df['main'].ffill()
			df.loc[~maincat,'BALANCE SHEET'] = df['main'] + ': ' + df['BALANCE SHEET']

			df['BALANCE SHEET'] = df['side'] + ':' + df['BALANCE SHEET'].apply(lambda x: re.sub(r'(-|,)', '', str(x).lower()))
			df['BALANCE SHEET'] = df['BALANCE SHEET'].apply(lambda x: re.sub(r'  ', ' ', x))

			df = df.set_index(df['BALANCE SHEET'])[['Unnamed: 4']]
			df = df.transpose()
			df.index = [pd.to_datetime(datetime.date(year, 1 + 3*(quarter-1), 1))]
			df['name'] = firms + ' Dealers'
			all_df = all_df.append(df)



for firms in ['Top 10', 'Top 11-25']:
	df = pd.read_excel('../input/SIFMA early_years.xlsx', sheetname=firms)
	df.loc[[('TOTAL LIABILITIES' in str(d)) and ('&' not in str(d)) for d in df['BALANCE SHEET']], 'side'] = 'Liabilities'
	df.loc[['TOTAL ASSETS' in str(d) for d in df['BALANCE SHEET']], 'side'] = 'Assets'
	df['side'] = df['side'].bfill()
	df = df.dropna(subset = [' ', 'BALANCE SHEET', 'side'])
	maincat = df['BALANCE SHEET'].apply(lambda x: not bool(re.match(r' *-', x)))
	df.loc[maincat, 'main'] = df['BALANCE SHEET']
	df['main'] = df['main'].ffill()
	df.loc[~maincat,'BALANCE SHEET'] = df['main'] + ': ' + df['BALANCE SHEET']

	df['BALANCE SHEET'] = df['side'] + ':' + df['BALANCE SHEET'].apply(lambda x: re.sub(r'(-|,)', '', str(x).lower()))
	df['BALANCE SHEET'] = df['BALANCE SHEET'].apply(lambda x: re.sub(r'  ', ' ', x))
	df.columns = df.iloc[0]
	df = df.iloc[1:]
	df = df.set_index(df.columns[0])
	df = df.transpose()
	df = df.iloc[:-2]
	df['name'] = firms + ' Dealers'
	df['qt_dt'] = [pd.to_datetime(datetime.date(2000 + int(r.split(':')[0]), 1 + 3*(int(r.split(':')[1][1])-1), 1)) for r in df.index]
	df = df.set_index('qt_dt')
	all_df = all_df.append(df)

beta_formula_in = ['Liabilities:bank loans payable', 'Liabilities:obligation to return securities', \
	'Liabilities:payables brokers dealers clearing organizations', 'Liabilities:repurchase agreements', 'Liabilities:securities sold short',\
	'Liabilities:subordinated liabilities: securities borrowed']
beta_formula_unc = []

liabs_out = ['Liabilities:payables to customers', 'Liabilities:payables to noncustomers', \
'Liabilities:accounts payable & accrued liabilities', 'Liabilities:notes & mortgages payable', \
'Liabilities:subordinated liabilities', 'Liabilities:special liabilities']


assets_in_unc = ['Assets:securities & spot commodities owned: bankers acceptances cds commercial paper', 'Assets:securities & spot commodities owned: corporate obligations',\
	 'Assets:securities & spot commodities owned: stocks & warrants', 'Assets:securities & spot commodities owned: other securities', \
	 'Assets:securities & spot commodities owned: debt securities (focus iia only)', 'Assets:securities owned not readily marketable',\
	  'Assets:other investments not readily marketable', 'Assets:other assets']
assets_in_formula = ['Assets:receivables brokers dealers clearing organizations', 'Assets:reverse repurchase agreements', 'Assets:securities & spot commodities owned: options',\
	 'Assets:securities & spot commodities owned: arbitrage']



breakdown = all_df[(all_df.index == pd.to_datetime('2016-10-01'))].sum().to_frame().transpose()
breakdown = breakdown[beta_formula_in + liabs_out + ['Liabilities:total liabilities']]
for v in beta_formula_in + liabs_out + ['Liabilities:total liabilities']:
		breakdown[v] = 100 * breakdown[v]/breakdown['Liabilities:total liabilities']

breakdown_renames = {'Liabilities:bank loans payable': 'Bank Loans Payable', \
 'Liabilities:obligation to return securities': 'Obligation to Return Securities', \
	'Liabilities:payables brokers dealers clearing organizations': 'Payables to BDs, Clearing' , \
 'Liabilities:repurchase agreements': 'Repurchase Agreements' , \
 'Liabilities:securities sold short': 'Securities Sold Short' , \
	'Liabilities:subordinated liabilities: securities borrowed': 'Securities Borrowed' , \
 'Liabilities:payables to customers': 'Payables to Customers' , \
 'Liabilities:payables to noncustomers': 'Payables to Non-Customers' , \
'Liabilities:accounts payable & accrued liabilities': 'Accounts Payable and Accrued Liabilities' , \
 'Liabilities:notes & mortgages payable': 'Notes and Mortgages' , \
'Liabilities:subordinated liabilities': 'Subordinated Liabilities' , \
 'Liabilities:special liabilities': 'Special Liabilities'}
breakdown = breakdown.drop('Liabilities:total liabilities', axis=1)
breakdown = breakdown.rename(columns = breakdown_renames).transpose()
breakdown = breakdown.astype(float).round(decimals=1)
#breakdown.columns = ['value']
#breakdown = breakdown.reset_index()
#breakdown = breakdown['index', 'value']
pd.options.display.max_colwidth = 100

f = open('../output/asset_breakdown_dealertop25.tex', 'wb')
f.write(breakdown.to_latex())
f.close()

all_df = all_df.sort_index()
all_df = all_df.rename(columns={'Assets:total assets': 'BHCK2170', 'Liabilities:total liabilities': 'BHCK2948'})
all_df[['BHCK2170', 'BHCK2948']] = 1000*all_df[['BHCK2170', 'BHCK2948']]


all_df['liab_in'] = 1000*all_df[beta_formula_in].sum(axis=1)
all_df['liab_in_unc'] = 0
all_df['liabs_in_frac'] = 1000*all_df[beta_formula_in].sum(axis=1)/all_df['BHCK2948']
all_df['asset_in'] = 1000*all_df[assets_in_formula].sum(axis=1)
all_df['asset_in_unc'] = 1000*all_df[assets_in_unc].sum(axis=1)

all_df = all_df[['BHCK2170', 'BHCK2948', 'name', 'liabs_in_frac', 'liab_in', 'liab_in_unc', 'asset_in', 'asset_in_unc']]
all_df.to_csv('../input/all.csv')
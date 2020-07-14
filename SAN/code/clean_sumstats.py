import re

renames = {'JPMORGAN': 'JP Morgan', 'AMER': 'America', 'FNCL': 'Financial', 'SVC': 'Services', 'NY': 'NY', 'FC': 'Finance Company', 'BK': 'Bank', 'BC': 'Bank', 'GRP': 'Group', 'OF': 'of', 'BBT': 'BBT'}

for j in ['asset', 'liab']:
	if j == 'asset':
		tag = 'Assets'
	elif j == 'liab':
		tag == 'Liabilities'

	txt = {'in': open('../output/' + j + '_in.tex', 'rb').read() , 'out': open('../output/' + j + '_out.tex', 'rb').read()}
	txt['in'] = re.sub(r'\s+&\s+\\\\', r'', txt['in'])
	txt['out'] = re.sub(r'\s+&\s+\\\\', r'', txt['out'])
	max_lines = max([len(re.findall(r'\\\\', txt['in'])), len(re.findall(r'\\\\', txt['out']))])
	for side in txt.keys():
		curr = txt[side]
		curr = re.sub(r'(\\begin{table.*\\begin{tabular)', r'\\begin{tabular', curr, flags=re.DOTALL)
		curr = re.sub(r'%', r'\\%', curr)
		curr = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', curr, flags=re.DOTALL)
		to_add = (max_lines - len(re.findall(r'\\\\', curr)))
		to_add = '\n\\\\addlinespace\n&\\\\\\' * to_add
		curr = re.sub(r'\\addlinespace.*?\n.*?\\%', to_add + '\n' +  r'\\addlinespace\n\%', curr)
		curr = re.sub(r'\\begin\{tabular\}', r'\\begin{tabular*}{.95\\textwidth}', curr)
		curr = re.sub(r'\\end\{tabular', r'\\end{tabular*', curr)
		curr = re.sub(r'&\\multicolumn\{1', r'\\multicolumn{2', curr)
		curr = re.sub(r'BHC(.*?)Financial(.*?)\}', r'\\begin{tabular}{@{}c@{}}BHC\1\\\\Financial\2\\end{tabular}} ', curr)
		curr = re.sub(r'\\bottomrule', r'', curr)
		curr = re.sub(r'\n\\% of BHC', r'\n\\bottomrule\n\% of BHC', curr)
		f = open('../output/' + j + '_' + side + '_formatted.tex', 'wb')
		f.write(curr)
		f.close()

	#txt_in = open('../output/' + j + '_in.tex', 'rb').read()
	#txt_out = open('../output/' + j + '_out.tex', 'rb').read()
	#max_lines = max([len(re.findall(r'\\\\', txt_out)), len(re.findall(r'\\\\', txt_in))])
	#txt_in = re.sub(r'(\\begin{table.*\\begin{tabular)', r'\\begin{tabular', txt_in, flags=re.DOTALL)
	#txt_out = re.sub(r'(\\begin{table.*\\begin{tabular)', r'\\begin{tabular', txt_out, flags=re.DOTALL)

	#txt_in = re.sub(r'%', r'\\%', txt_in)
	#txt_out = re.sub(r'%', r'\\%', txt_out)

	#txt_in = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', txt_in, flags=re.DOTALL)
	#txt_out = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', txt_out, flags=re.DOTALL)

	#to_add = (max_lines - len(re.findall(r'\\\\', txt_out)))/2
	#to_add = '\n\\\\addlinespace\n\n&\\\\\\\n&\\\\\\' * to_add
	#txt_out = re.sub(r'\\addlinespace.*?\n.*?\\%', to_add + '\n' +  r'\\addlinespace\n\%', txt_out)
"""
	to_add = (max_lines - len(re.findall(r'\\\\', txt_in)))/2
	to_add = '\n\\\\addlinespace\n\n&\\\\\\\n&\\\\\\' * to_add
	txt_in = re.sub(r'\\addlinespace.*?\n\\%', to_add + '\n' +  r'\\addlinespace\n\%', txt_in)

	txt_in = re.sub(r'\\begin\{tabular\}', r'\\begin{tabular*}{.95\\textwidth}', txt_in)
	txt_out = re.sub(r'\\begin\{tabular\}', r'\\begin{tabular*}{.95\\textwidth}', txt_out)

	txt_in = re.sub(r'\\end\{tabular', r'\\end{tabular*', txt_in)
	txt_out = re.sub(r'\\end\{tabular', r'\\end{tabular*', txt_out)

	txt_in = re.sub(r'&\\multicolumn\{1', r'\\multicolumn{2', txt_in)
	txt_out = re.sub(r'&\\multicolumn\{1', r'\\multicolumn{2', txt_out)

	txt_in = re.sub(r'BHC(.*?)Financial(.*?)\}', r'\\begin{tabular}{@{}c@{}}BHC\1\\\\Financial\2\\end{tabular}} ', txt_in)
	txt_out = re.sub(r'BHC(.*?)Financial(.*?)\}', r'\\begin{tabular}{@{}c@{}}BHC\1\\\\Financial\2\\end{tabular}} ', txt_out)

	txt_in = re.sub(r'\\bottomrule', r'', txt_in)
	txt_out = re.sub(r'\\bottomrule', r'', txt_out)

	txt_in = re.sub(r'\n\\% of BHC', r'\n\\bottomrule\n\% of BHC', txt_in)
	txt_out = re.sub(r'\n\\% of BHC', r'\n\\bottomrule\n\% of BHC', txt_out)

	f = open('../output/' + j + '_in_formatted.tex', 'wb')
	f.write(txt_in)
	f.close()


	f= open('../output/' + j + '_out_formatted.tex', 'wb')
	f.write(txt_out)
	f.close()
"""

#Compustat samples stuff
headers = {'dealers' : 'Security Broker and Dealer', 'insurance' : 'Insurance Company', 'reit': 'Real Estate Investment Trust', 'other': 'Other Financial Firm', 'dealers_top10': 'Top 10 Dealers by Assets', 'dealers_top25': 'Top 11-25 Dealers by Assets'}
for tag in ['insurance', 'reit', 'other', 'dealers_top10', 'dealers_top25']:
	curr = open('../output/' + tag + '_sample.tex', 'rb').read()
	curr = re.sub(r'\\begin{table.*\\begin{tabular', r'\\begin{tabular', curr, flags=re.DOTALL)
	curr = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', curr, flags=re.DOTALL)
	curr = re.sub(r'\\bottomrule', r'', curr)
	curr = re.sub(r'\nNumber of Firms(.*?)\.00', r'\n\\bottomrule\nNumber of Firms\1', curr)
	curr = re.sub(r'\n\s+&\\multicolumn', r'\n' + headers[tag] + r'&\multicolumn' , curr)
	curr = re.sub(r'-', r'', curr)
	curr = re.sub(r'([A-Z]+)', lambda m: m.group(1).title() if m.group(1) not in renames.keys() else renames[m.group(1)], curr)
	f = open('../output/' + tag + '.tex', 'wb')
	f.write(curr)
	f.close()

"""
dealers = open('../output/dealers_sample.tex', 'rb').read()
insurance = open('../output/insurance_sample.tex', 'rb').read()
reits = open('../output/reit_sample.tex', 'rb').read()

dealers = re.sub(r'\\begin{table.*\\begin{tabular', r'\\begin{tabular', dealers, flags=re.DOTALL)
insurance = re.sub(r'\\begin{table.*\\begin{tabular', r'\\begin{tabular', insurance, flags=re.DOTALL)
reits = re.sub(r'\\begin{table.*\\begin{tabular', r'\\begin{tabular', reits, flags=re.DOTALL)
dealers = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', dealers, flags=re.DOTALL)
insurance = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', insurance, flags=re.DOTALL)
reits = re.sub(r'\\end{tabular}.*\\end{table}.*', r'\\end{tabular}', reits, flags=re.DOTALL)


dealers = re.sub(r'\\bottomrule', r'', dealers)
insurance = re.sub(r'\\bottomrule', r'', insurance)
reits = re.sub(r'\\bottomrule', r'', reits)

dealers = re.sub(r'\nNumber of Firms(.*?)\.00', r'\n\\bottomrule\nNumber of Firms\1', dealers)
insurance = re.sub(r'\nNumber of Firms(.*?)\.00', r'\n\\bottomrule\nNumber of Firms\1', insurance)
reits = re.sub(r'\nNumber of Firms(.*?)\.00', r'\n\\bottomrule\nNumber of Firms\1', reits)
f = open('../output/dealers.tex', 'wb')
f.write(dealers)
f.close()


f= open('../output/insurance.tex', 'wb')
f.write(insurance)
f.close()

f= open('../output/reit.tex', 'wb')
f.write(reits)
f.close()
"""
keys = {'insurance': 'INSURANCE', 'dealer': 'DEALERS', 'reit': 'REIT', 'credit_union': 'CREDIT UNIONS', 'finance': 'FINANCE COMPANIES', 'funding': 'FUNDING CORPORATIONS', 'abs': 'ABS ISSUERS', 'dealertop25': 'DEALERS25'}
orders = {k: {} for k in keys.keys()}
#insurance_order = {}
#dealer_order = {}
#reit_order = {}
master_tables = open('../paper/master_table_assets_wffunds.tex', 'rb').read()

percents_text = {k: open('../output/asset_breakdown_' + k + '.tex', 'rb').read() for k in keys.keys()}

for k in keys.keys():
	categories = [s.strip() for s in sorted(re.findall(r'\} (.*?)&.*?--.*?' + keys[k] + '.*?--', master_tables), key=lambda x: -1 * len(x))]
	for cat in categories:
		num = re.search(r'' + cat + r'.*?&(.*?)\\', percents_text[k])
		if num:
			master_tables = re.sub(r'(\} ' + cat + r'.*?)--.*?' + keys[k] + '.*?--', r'\g<1>' + num.groups()[0].strip(), master_tables)
			orders[k][cat] = float(num.groups()[0].strip())
#perc_ins_ffunds = open('../output/asset_breakdown_insurance.tex', 'rb').read()
#perc_dealer_ffunds = open('../output/asset_breakdown_dealer.tex', 'rb').read()
#perc_reit_ffunds = open('../output/asset_breakdown_reit.tex', 'rb').read()



firm_specific = open('../output/large_instit_breakdowns.tex', 'rb').read()
firm_specific = re.sub(r'(.*\n)+(\\begin\{tabular(.*\n)+\\end\{tabular\}.*)(.*\n)+', r'\2', firm_specific)
firm_specific = re.sub(r'([A-Z]+)', lambda m: m.group(1).title() if m.group(1) not in renames.keys() else renames[m.group(1)] , firm_specific)
firm_specific = re.sub(r'.*\\multicolumn.*\n', r'& & \multirow{2}{*}{\shortstack[c]{Financial \\\\ Connectivity}} & &\\\\\n', firm_specific)
firm_specific = re.sub('&Financial Connectivity&', r'& &', firm_specific)
f = open('../output/large_instit_breakdowns_clean.tex', 'wb')
f.write(firm_specific)
f.close()

"""
categories_insurance = [s.strip() for s in sorted(re.findall(r'\} (.*?)&.*?--.*?INSURANCE.*?--', master_tables), key=lambda x: -1 * len(x))]
for cat in categories_insurance:
	num = re.search(r'' + cat + r'.*?&(.*?)\\', perc_ins_ffunds)
	if num:
		master_tables = re.sub(r'(\} ' + cat + r'.*?)--.*?INSURANCE.*?--', r'\g<1>' + num.groups()[0].strip(), master_tables)
		insurance_order[cat] = float(num.groups()[0].strip())

categories_dealers = [s.strip() for s in sorted(re.findall(r'\} (.*?)&.*?--.*?DEALERS.*?--', master_tables), key=lambda x: -1 * len(x))]
for cat in categories_dealers:
	num = re.search(r'' + cat + r'.*?&(.*?)\\', perc_dealer_ffunds)
	if num:
		master_tables = re.sub(r'(\} ' + cat + r'.*?)--.*?DEALERS.*?--', r'\g<1>' + num.groups()[0].strip(), master_tables)
		dealer_order[cat] = float(num.groups()[0].strip())

categories_reit = [s.strip() for s in sorted(re.findall(r'\} (.*?)&.*?--.*?REIT.*?--', master_tables), key=lambda x: -1 * len(x))]
for cat in categories_reit:
	num = re.search(r'' + cat + r'.*?&(.*?)\\', perc_reit_ffunds)
	if num:
		master_tables = re.sub(r'(\} ' + cat + r'.*?)--.*?REIT.*?--', r'\g<1>' + num.groups()[0].strip(), master_tables)
		reit_order[cat] = float(num.groups()[0].strip())
"""
f = open('../output/master_table_assets_wffunds_inserts.tex', 'wb')
f.write(master_tables)
f.close()

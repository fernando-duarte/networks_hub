#!/data/apps/Anaconda2-4.3.1/bin/python

import subprocess
import pexpect
import sys
import time
import pandas as pd

wrds_username = ''
passwd = ''
#Also change the contents of batch_pull.sh, as needed

#Must create .pgpass file in your WRDS cloud home directory, it you don't already have it (see https://wrds-www.wharton.upenn.edu/pages/support/programming-wrds/python-programming-wrds/)
#Make this directory in your WRDS cloud home directory, if you don't already have it
pull_directory = '/home/frb-ny/' + wrds_username + '/'

child = pexpect.spawn('scp WRDS_query.py ' + wrds_username + '@wrds-cloud.wharton.upenn.edu:' + pull_directory, timeout=30)
child.expect(wrds_username + '@wrds-cloud.wharton.upenn.edu\'s password:')
child.sendline(passwd)
child.expect(pexpect.EOF)
child.close()

child = pexpect.spawn('scp batch_pull.sh ' + wrds_username + '@wrds-cloud.wharton.upenn.edu:' + pull_directory, timeout=30)
child.expect(wrds_username + '@wrds-cloud.wharton.upenn.edu\'s password:')
child.sendline(passwd)
child.expect(pexpect.EOF)
child.close()

child = pexpect.spawn('ssh ' + wrds_username + '@wrds-cloud.wharton.upenn.edu qsub ' + pull_directory + 'batch_pull.sh', timeout=60)
child.expect(wrds_username + '@wrds-cloud.wharton.upenn.edu\'s password:')
child.sendline(passwd)
child.expect(pexpect.EOF)
child.close()

time.sleep(15)

child = pexpect.spawn('scp ' + wrds_username + '@wrds-cloud.wharton.upenn.edu:' + pull_directory + 'permco_cusip.csv ../temp/', timeout=30)
child.expect(wrds_username + '@wrds-cloud.wharton.upenn.edu\'s password:')
child.sendline(passwd)
child.expect(pexpect.EOF)
child.close()
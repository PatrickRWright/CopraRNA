import sys
IntaRNA_result = sys.argv[1]
#print IntaRNA_result

with open(IntaRNA_result) as file:
    IntaRNA_lines = file.readlines()

#print IntaRNA_lines[1]

backgroundList = []
# go through IntaRNA output line by line and extract Entrez Gene Ids
for i in range(1,len(IntaRNA_lines)): # range omits the right boundary
    # split
    curr_line = IntaRNA_lines[i]
    split = curr_line.split(";")
    entrezID = split[9]
    backgroundList.append(entrezID)

backgroundList = map(str,backgroundList)

print "background:" + str(len(backgroundList))

inputList = backgroundList[0:100]
print "input:" + str(len(inputList))

inputIds = ",".join(inputList)
#print inputIds
backgroundIds = ",".join(backgroundList)
#print backgroundIds

#sys.exit()

sys.path.append('../')

import logging
import traceback as tb
import suds.metrics as metrics
#from tests import *
from suds import *
from suds.client import Client
from datetime import datetime

errors = 0

#setup_logging()

#logging.getLogger('suds.client').setLevel(logging.DEBUG)

url = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
print 'url=%s' % url
# create a service client using the wsdl.
client = Client(url)
ws = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'
client.wsdl.services[0].setlocation(ws)

exit
#
# print the service (introspection)
#
print client

#authenticate user email 
print client.service.authenticate('patrickrw@gmx.net')

# add top 100 predicted
idType = 'ENTREZ_GENE_ID'
listName = 'make_up'
listType = 0
print client.service.addList(inputIds, idType, listName, listType)

# add background
#listNameBack = 'back'
#listTypeBack = 1
#print client.service.addList(backgroundIds, idType, listNameBack, listTypeBack)

print client.service.getDefaultCategoryNames()

#getTermClusterReport
thd = 1
ct = 1
print client.service.getChartReport(thd,ct)


#!/usr/bin/env python

import sys
CopraRNA_result = sys.argv[1]
# amount of top predictions to enrich with
enrichment_count = int(sys.argv[2])
#print CopraRNA_result

with open(CopraRNA_result) as file:
    CopraRNA_lines = file.readlines()

backgroundList = []
# go through CopraRNA output line by line and extract Entrez Gene Ids
for i in range(1,len(CopraRNA_lines)): # range omits the right boundary
    # split
    curr_line = CopraRNA_lines[i]
    split = curr_line.split(",")
    orgOfIntEntry = split[3]
    if orgOfIntEntry:
        splitOrgOfInt = orgOfIntEntry.split(":")
        entrezID = splitOrgOfInt[1][:-1]
        backgroundList.append(entrezID)

backgroundList = list(map(str,backgroundList))

print ("background:" + str(len(backgroundList)))

inputList = backgroundList[0:enrichment_count]
print ("input:" + str(len(inputList)))

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
print ('url=%s' % url)
# create a service client using the wsdl.
client = Client(url)
ws = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'
client.wsdl.services[0].setlocation(ws)

exit
#
# print the service (introspection)
#
print (client)

#authenticate user email 
print (client.service.authenticate('patrickrw@gmx.net'))

# add enrichment_count (amount) predicted
idType = 'ENTREZ_GENE_ID'
listName = 'make_up'
listType = 0
print (client.service.addList(inputIds, idType, listName, listType))

# add background
listNameBack = 'back'
listTypeBack = 1
print (client.service.addList(backgroundIds, idType, listNameBack, listTypeBack))

print (client.service.getDefaultCategoryNames())

##getChartReport
#thd=0.1
#count = 2
#print client.service.getChartReport(thd, count)

#getTermClusterReport
overlap = 3
initialSeed = 2
finalSeed = 2
linkage = 0.5
kappa = 85
print (client.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa))



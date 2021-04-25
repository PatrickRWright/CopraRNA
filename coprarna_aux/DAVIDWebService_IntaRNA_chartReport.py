#!/usr/bin/env python

import sys

if len(sys.argv) < 3:
   print("\ncall: DAVIDWebService_IntaRNA_chartReport.py <CSV-file-with-entrex_gene_id-col> <enrich_count>\n");
   exit(-1);

IntaRNA_result = sys.argv[1]
enrich_count = int(sys.argv[2])
#print IntaRNA_result

with open(IntaRNA_result) as file:
    IntaRNA_lines = file.readlines()

# identify entrez_gene_id column index
entrezIdCol = IntaRNA_lines[0].split(";").index("entrez_gene_id")

#print IntaRNA_lines[1]

backgroundList = []
# go through IntaRNA output line by line and extract Entrez Gene Ids
for i in range(1,len(IntaRNA_lines)): # range omits the right boundary
    # extract entrezID for this line
    entrezID = IntaRNA_lines[i].split(";")[entrezIdCol]
    # store id
    backgroundList.append(entrezID)

backgroundList = list(map(str,backgroundList))

print ("background:" + str(len(backgroundList)))

inputList = backgroundList[0:enrich_count]
print ("input:" + str(len(inputList)))

inputIds = ",".join(inputList)
#print inputIds
backgroundIds = ",".join(backgroundList)
#print backgroundIds

#sys.exit()

#sys.path.append('../')

#import logging
import traceback as tb
import suds.metrics as metrics
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

#exit
#
# print the service (introspection)
#
print (client)

#authenticate user email 
print (client.service.authenticate('rna@informatik.uni-freiburg.de'))

print(str(inputIds))

# add enrich_count (amount) predicted
idType = 'ENTREZ_GENE_ID'
listName = 'make_up'
listType = 0
print (client.service.addList(inputIds, idType, listName, listType))

print (client.service.getDefaultCategoryNames())

thd = 1
ct = 1
print (client.service.getChartReport(thd,ct))


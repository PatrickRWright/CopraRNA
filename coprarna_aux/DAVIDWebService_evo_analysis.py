import sys
evo_prep_result = sys.argv[1]

with open(evo_prep_result) as file:
    evo_prep_lines = file.readlines()

# hash for each RefSeq ID pointing at the concatenated list of GIDs
header = evo_prep_lines[0]
header = header.rstrip('\n') 
split_header = header.split(",")

refid_to_gids = {}

for i in range(1,len(split_header)):
    refid_to_gids[i] = split_header[i];

# go through evo prep output line by line and extract Entrez Gene Ids
for i in range(1,151): # top 150 lines -> in this case that does not mean we retrieve 150 GIDs
    # split
    curr_line = evo_prep_lines[i]
    curr_line = curr_line.rstrip('\n')
    split = curr_line.split(",")
    for j in range(1,len(split)):
        refid_to_gids[j] = refid_to_gids[j] + "," + split[j]    

#print refid_to_gids[1]
#print refid_to_gids[2]
#print refid_to_gids[3]
#print refid_to_gids[4]
#print len(refid_to_gids)

for i in range(1,(len(refid_to_gids)+1)):
    
    inputIds_split = refid_to_gids[i].split(",")
    RefSeqID = inputIds_split.pop(0)
    print RefSeqID
    inputIds_split = filter(None, inputIds_split)
    print inputIds_split
    print len(inputIds_split)
    inputIds = ",".join(inputIds_split)
    #sys.exit()

    sys.path.append('../')

    import logging
    import traceback as tb
    import suds.metrics as metrics
    from suds import *
    from suds.client import Client
    from datetime import datetime

    errors = 0

    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
    print 'url=%s' % url
    # create a service client using the wsdl.
    client = Client(url)
    ws = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/'
    client.wsdl.services[0].setlocation(ws)

    exit
    #
    # print the service (introspection)
    #
    print client

    #authenticate user email 
    print client.service.authenticate('patrickrw@gmx.net')

    # add input
    idType = 'ENTREZ_GENE_ID'
    listName = 'make_up'
    listType = 0
    print client.service.addList(inputIds, idType, listName, listType)

    print client.service.getDefaultCategoryNames()
    write_file = RefSeqID + '-chartReport.txt';
    f = open(write_file, 'w')
    #getTermClusterReport
    thd = 1
    ct = 1
    chart_report = str(client.service.getChartReport(thd,ct))
    f.write(chart_report)
    f.close()



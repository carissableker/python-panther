#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import os
import sys
import requests
from bs4 import BeautifulSoup
import argparse



# ### Annotation Dataset: 
# 
# | Option           | Annotation Dataset        |
# |------------------|--------------------------|
# | `pathway`        | PANTHER Pathways</option> |
# | `panther_mf`     | PANTHER GO-Slim Molecular Function|
# | `panther_bp`     | PANTHER GO-Slim Biological Process|
# | `panther_cc`     | PANTHER GO-Slim Cellular Component|
# | `panther_pc`     | PANTHER Protein Class|
# | `fullgo_mf_comp` | GO molecular function complete|
# | `fullgo_bp_comp` | GO biological process complete|                                  
# | `fullgo_cc_comp` | GO cellular component complete|
# | `reactome`       | Reactome pathways
#                     

def panther_api_overrepresentation(inputfile, organism, annotation_option, test_type):
    url = 'http://pantherdb.org/webservices/bdds/overrep.jsp?'

    files = {'organism':       (None, organism), 
             'input':          (None, open(inputfile, 'r')),
             'enrichmentType': (None, annotation_option), 
             'test_type':      (None, test_type)}

    response = requests.post(url=url, files=files)
    soup = BeautifulSoup(response.text, "html5lib")
    response.close()

    max_sessions_error = "Maximum number of sessions exceeded, please try later"
    if soup.text == max_sessions_error:
        print(max_sessions_error)
        sys.exit()

    n = soup.find_all(name='a', attrs={'href':'/tools/gxIdsList.do?list=upload_1&organism=%s'%organism})[0].text
    not_n = soup.find_all(name='a', attrs={'href':'/tools/unmappedBinom.jsp?listName=upload_1'})[0].text
    N = soup.find_all(name='a', attrs={'href':'/tools/gxIdsList.do?reflist=1'})[0].text

    print('Reference size:        ', N)
    print('Number IDs mapped:     ', n) 
    print('Number IDs not mapped: ', not_n) 

    if 'No statistically significant results.' in soup.text:
        print('No statistically significant results.')
        return 

    GO_terms = {}
    tag_list_term = soup.find_all(name='td', attrs={'class':'fixHeader', 'nowrap':''})
    tag_list_result = soup.find_all(name='td', attrs={'class':'tools'}) 
    i = 0
    for tag in tag_list_term:
        D = {}
        for x in tag.children:
            if x.name == 'a':
                GO_name = x.text
                GO_id = x.attrs['href'].split('/')[-1]
                break
        
        # numbers
        go_href = "/tools/gxIdsList.do?acc=%s&list=upload_1&organism=%s"%(GO_id, organism)    
        D['# in list'] = int(soup.find_all('a', attrs={'href':go_href})[0].text.strip())
        ref_href = '/tools/gxIdsList.do?acc=%s&reflist=1'%GO_id
        D['# in reference'] = int(soup.find_all('a', attrs={'href':ref_href})[0].text.strip())

        enrichment_result = tag_list_result[i:i+5]

        D['# expected in list']  = enrichment_result[0].text
        D['fold_enrichment']     = enrichment_result[1].text
        D['direction']           = enrichment_result[2].text
        D['pvalue']              = enrichment_result[3].text
        D['FDR']                 = enrichment_result[4].text

        D['name'] = GO_name
        GO_terms[GO_id] = D
        
        i += 5

    df = pd.DataFrame(GO_terms).T
    df = df[['name', '# in reference', '# in list', '# expected in list', 'fold_enrichment', 'direction', 'pvalue', 'FDR']]
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PANTHER overenrichment test of a gene list.')
    parser.add_argument("inputfile", type=str, help="Gene list file. One ID per line.")
    parser.add_argument("outputfile", type=str, help="File to save results")
    parser.add_argument("--organism", type=str, help="Organism for reference/background", default='Homo sapiens')
    parser.add_argument("--test_type", type=str, help="One of FISHER or BINOMIAL", default='FISHER')
    parser.add_argument("--annotation_option", type=str, help="Annotation option, see table in code", default='fullgo_bp_comp')
    args = parser.parse_args()

    result = panther_api_overrepresentation(args.inputfile, args.organism, args.annotation_option, args.test_type)
    if result is not None:
        result.to_csv(args.outputfile, sep='\t')



# expected = # number in list*(# number in reference/reference size)

#!/bin/python
#!/usr/bin/env python3
# Originally written by Keith Jolley and adapted by ycode-sh for PUBMLST profiling 

import sys, requests, argparse, base64
parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', type=str, default='contigs.fasta', help='assembly contig filename (FASTA format)')
args = parser.parse_args()

def main():
    uri = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence'
    with open(args.file, 'r') as x: 
        fasta = x.read()
    payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}'
    response = requests.post(uri, data=payload)
    if response.status_code == requests.codes.ok:
        data = response.json()
        try: 
            data['taxon_prediction']
        except KeyError:
            print("No match")
            sys.exit(0)
        for match in data['taxon_prediction']:
                print("Rank: " + match['rank'])
                print("Taxon:" + match['taxon'])
                print("Support:" + str(match['support']) + "%")
                print("Taxonomy:" + match['taxonomy'] + "\n")
            
    else:
        print(response.text)

if __name__ == "__main__":
    main()

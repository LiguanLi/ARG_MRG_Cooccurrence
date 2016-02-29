#! /usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import getopt

#file_name1=sys.argv[1]
file_name1=raw_input("Enter the full name of genebank file: ")

gb = SeqIO.parse(file_name1, "genbank")
output_handle1=open(file_name1+".ffn","w")
output_handle2=open(file_name1+".faa","w")
output_handle3=open(file_name1+".info.tab","w")

for record in gb:
    print "Processing GenBank record %s" % record.id
    for feature in record.features :
        if(feature.type == "source"):
            taxaid = feature.qualifiers['db_xref']
            if len(taxaid) > 1 :
                output_handle3.write("%s\t%s\t%s\t\n" % (record.id, record.description, taxaid[1]))
            else:
                output_handle3.write("%s\t%s\t%s\t\n" % (record.id, record.description, taxaid[0]))
    for feature in record.features :
        if(feature.type == "CDS"):
            try:
                ID = feature.qualifiers['db_xref'][0]
                desc = feature.qualifiers['protein_id'][0]
                product = feature.qualifiers['product'][0]
                locus = feature.qualifiers['locus_tag'][0]
                nt_seq = feature.extract(record.seq)
                strand=feature.strand
                start=location.start.position
                end=location.end.position
                try:
                        assert len(feature.qualifiers['translation'])==1
                        aa_seq = feature.qualifiers['translation'][0]
                except AssertionError:
                        print ID,'no amni acids found!'
                        aa_seq=''  
                output_handle1.write(">%s %s %s %s %s %s %s\n%s\n" % (ID, desc, product, locus, strand, start, end, nt_seq))
                output_handle2.write(">%s %s %s %s %s %s %s\n%s\n" % (ID, desc, product, locus, strand, start, end, aa_seq))
            except:
                continue

print 'Retrieving whole genome sequences!'
output_handle1.close()
output_handle2.close()
output_handle3.close()
print 'Done!'

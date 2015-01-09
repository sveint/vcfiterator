# vcfiterator
VCF iterator, Python module for iterating over a VCF file, parsing and generating data structures.

This module is still very much a WIP, but it should work.

## Purpose

Parsing VCF files into a more machine usable format can be difficult, especially when dealing with INFO fields from annotation software.

The main purpose for this module is to be able to parse VCF files and generate usable JSON from it, with no further parsing necessary.
It can also be useful for any task where the annotation data needs to be readily available.

## Features

It currently supports output from the following annotation software:

1. VEP
2. snpEff

You can easily add support for more software, pull requests are welcome.

It does not yet support tabix indexing or compressed VCF files.

## Usage

```
v = VcfIterator(path)

# Print all VEP info for all variants
for variant in v.iter():
    alleles = variant['ALT']
    for allele in alleles:
        print allele
        print variant['INFO'][allele]['CSQ']
```

## Example output

In short, it takes a VCF line like the following:
```
5   179390472   .   C   A   21.77   LowQual AC=2;AF=1.0;AN=2;DP=2;FS=0.0;MLEAC=2;MLEAF=1.0;MQ=60.0;MQ0=0;QD=10.88;CSQ=A|55819|NM_018434.5|Transcript|stop_gained&splice_region_variant|1659|1243|415|E/*|Gag/Tag||1||-1|protein_coding|YES||NP_060904.2|||8/9|||NM_018434.5:c.1243G>T|NP_060904.2:p.Glu415Ter||||||||||||||,A|ENSG00000113269|ENST00000519708|Transcript|downstream_gene_variant|||||||1|3511|-1|retained_intron||||||||||||||||||||||||,A|55819|NM_001280801.1|Transcript|intron_variant|||||||1||-1|protein_coding|||NP_001267730.1||||7/7||NM_001280801.1:c.1150+3334G>T|||||||||||||||,A|ENSG00000113269|ENST00000521389|Transcript|stop_gained&splice_region_variant|1659|1243|415|E/*|Gag/Tag||1||-1|protein_coding|YES|CCDS4451.1|ENSP00000430237|||8/9|||ENST00000521389.1:c.1243G>T|ENSP00000430237.1:p.Glu415Ter||||||||||||||,A|ENSG00000113269|ENST00000522208|Transcript|intron_variant|||||||1||-1|protein_coding|||ENSP00000429509||||7/7||ENST00000522208.2:c.1150+3334G>T|||||||||||||||,A|ENSG00000113269|ENST00000521901|Transcript|splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant|472||||||1||-1|retained_intron||||||1/2|||ENST00000521901.1:n.472G>T|||||||||||||||,A|ENSG00000249412|ENST00000510240|Transcript|upstream_gene_variant|||||||1|143|1|antisense|YES|||||||||||||||||||||||,A|ENSG00000113269|ENST00000520911|Transcript|splice_region_variant&3_prime_UTR_variant&NMD_transcript_variant|1599||||||1||-1|nonsense_mediated_decay|||ENSP00000430999|||8/9|||ENST00000520911.1:c.*762G>T|||||||||||||||,A|ENSG00000113269|ENST00000261947|Transcript|intron_variant|||||||1||-1|protein_coding||CCDS64340.1|ENSP00000261947||||7/7||ENST00000261947.4:c.1150+3334G>T|||||||||||||||;EFF=stop_gained(HIGH|NONSENSE|Gag/Tag|p.Glu415*/c.1243G>T|419|RNF130|protein_coding|CODING|ENST00000521389|8|1),splice_region_variant(LOW|||c.1243G>T|419|RNF130|protein_coding|CODING|ENST00000521389|8|1),splice_region_variant(LOW|||||RNF130|nonsense_mediated_decay|CODING|ENST00000520911|8|1),splice_region_variant(LOW|||n.472G>T||RNF130|retained_intron|CODING|ENST00000521901|1|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|7|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|4|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|3|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|9|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|6|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|5|1),sequence_feature[topological_domain:Cytoplasmic](LOW|||c.1243C>A|419|RNF130|protein_coding|CODING|ENST00000521389|8|1),3_prime_UTR_variant(MODIFIER||49801|n.*762G>T||RNF130|nonsense_mediated_decay|CODING|ENST00000520911|8|1),upstream_gene_variant(MODIFIER||143|||CTC-563A5.2|antisense|NON_CODING|ENST00000510240||1),downstream_gene_variant(MODIFIER||3511|||RNF130|retained_intron|CODING|ENST00000519708||1),intron_variant(MODIFIER|||c.1150+3334G>T|419|RNF130|protein_coding|CODING|ENST00000522208|7|1),intron_variant(MODIFIER|||c.1150+3334G>T|384|RNF130|protein_coding|CODING|ENST00000261947|7|1),non_coding_exon_variant(MODIFIER|||n.472G>T||RNF130|retained_intron|CODING|ENST00000521901|1|1) GT:AD:DP:GQ:PL  1/1:0,2:2:6:49,6,0
```
and turns it into a data structure:

```
{'ALT': ['A'],
 'CHROM': '5',
 'FILTER': 'LowQual',
 'ID': '.',
 'INFO': {'A': {'AC': '2',
                'AF': '1.0',
                'CSQ': [{'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'Amino_acids': 'E/*',
                         'BIOTYPE': 'protein_coding',
                         'CANONICAL': 'YES',
                         'CDS_position': '1243',
                         'Codons': 'Gag/Tag',
                         'Consequence': 'stop_gained&splice_region_variant',
                         'ENSP': 'NP_060904.2',
                         'EXON': '8/9',
                         'Feature': 'NM_018434.5',
                         'Feature_type': 'Transcript',
                         'Gene': '55819',
                         'HGVSc': 'NM_018434.5:c.1243G>T',
                         'HGVSp': 'NP_060904.2:p.Glu415Ter',
                         'Protein_position': '415',
                         'STRAND': -1,
                         'cDNA_position': '1659'},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'retained_intron',
                         'Consequence': 'downstream_gene_variant',
                         'DISTANCE': 3511,
                         'Feature': 'ENST00000519708',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'STRAND': -1},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'protein_coding',
                         'Consequence': 'intron_variant',
                         'ENSP': 'NP_001267730.1',
                         'Feature': 'NM_001280801.1',
                         'Feature_type': 'Transcript',
                         'Gene': '55819',
                         'HGVSc': 'NM_001280801.1:c.1150+3334G>T',
                         'INTRON': '7/7',
                         'STRAND': -1},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'Amino_acids': 'E/*',
                         'BIOTYPE': 'protein_coding',
                         'CANONICAL': 'YES',
                         'CCDS': 'CCDS4451.1',
                         'CDS_position': '1243',
                         'Codons': 'Gag/Tag',
                         'Consequence': 'stop_gained&splice_region_variant',
                         'ENSP': 'ENSP00000430237',
                         'EXON': '8/9',
                         'Feature': 'ENST00000521389',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'HGVSc': 'ENST00000521389.1:c.1243G>T',
                         'HGVSp': 'ENSP00000430237.1:p.Glu415Ter',
                         'Protein_position': '415',
                         'STRAND': -1,
                         'cDNA_position': '1659'},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'protein_coding',
                         'Consequence': 'intron_variant',
                         'ENSP': 'ENSP00000429509',
                         'Feature': 'ENST00000522208',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'HGVSc': 'ENST00000522208.2:c.1150+3334G>T',
                         'INTRON': '7/7',
                         'STRAND': -1},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'retained_intron',
                         'Consequence': 'splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant',
                         'EXON': '1/2',
                         'Feature': 'ENST00000521901',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'HGVSc': 'ENST00000521901.1:n.472G>T',
                         'STRAND': -1,
                         'cDNA_position': '472'},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'antisense',
                         'CANONICAL': 'YES',
                         'Consequence': 'upstream_gene_variant',
                         'DISTANCE': 143,
                         'Feature': 'ENST00000510240',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000249412',
                         'STRAND': 1},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'nonsense_mediated_decay',
                         'Consequence': 'splice_region_variant&3_prime_UTR_variant&NMD_transcript_variant',
                         'ENSP': 'ENSP00000430999',
                         'EXON': '8/9',
                         'Feature': 'ENST00000520911',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'HGVSc': 'ENST00000520911.1:c.*762G>T',
                         'STRAND': -1,
                         'cDNA_position': '1599'},
                        {'ALLELE_NUM': 1,
                         'Allele': 'A',
                         'BIOTYPE': 'protein_coding',
                         'CCDS': 'CCDS64340.1',
                         'Consequence': 'intron_variant',
                         'ENSP': 'ENSP00000261947',
                         'Feature': 'ENST00000261947',
                         'Feature_type': 'Transcript',
                         'Gene': 'ENSG00000113269',
                         'HGVSc': 'ENST00000261947.4:c.1150+3334G>T',
                         'INTRON': '7/7',
                         'STRAND': -1}],
                'EFF': [{'Amino_Acid_Change': 'p.Glu415*/c.1243G>T',
                         'Amino_Acid_length': 419,
                         'Codon_Change': 'Gag/Tag',
                         'Effect': 'stop_gained',
                         'Effect_Impact': 'HIGH',
                         'Exon_Rank': 8,
                         'Functional_Class': 'NONSENSE',
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243G>T',
                         'Amino_Acid_length': 419,
                         'Effect': 'splice_region_variant',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 8,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Effect': 'splice_region_variant',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 8,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'nonsense_mediated_decay',
                         'Transcript_ID': 'ENST00000520911'},
                        {'Amino_Acid_Change': 'n.472G>T',
                         'Effect': 'splice_region_variant',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 1,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'retained_intron',
                         'Transcript_ID': 'ENST00000521901'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 7,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 4,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 3,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 9,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 6,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 5,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'c.1243C>A',
                         'Amino_Acid_length': 419,
                         'Effect': 'sequence_feature[topological_domain:Cytoplasmic]',
                         'Effect_Impact': 'LOW',
                         'Exon_Rank': 8,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000521389'},
                        {'Amino_Acid_Change': 'n.*762G>T',
                         'Codon_Change': '49801',
                         'Effect': '3_prime_UTR_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Exon_Rank': 8,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'nonsense_mediated_decay',
                         'Transcript_ID': 'ENST00000520911'},
                        {'Codon_Change': '143',
                         'Effect': 'upstream_gene_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Gene_Coding': 'NON_CODING',
                         'Gene_Name': 'CTC-563A5.2',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'antisense',
                         'Transcript_ID': 'ENST00000510240'},
                        {'Codon_Change': '3511',
                         'Effect': 'downstream_gene_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'retained_intron',
                         'Transcript_ID': 'ENST00000519708'},
                        {'Amino_Acid_Change': 'c.1150+3334G>T',
                         'Amino_Acid_length': 419,
                         'Effect': 'intron_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Exon_Rank': 7,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000522208'},
                        {'Amino_Acid_Change': 'c.1150+3334G>T',
                         'Amino_Acid_length': 384,
                         'Effect': 'intron_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Exon_Rank': 7,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'protein_coding',
                         'Transcript_ID': 'ENST00000261947'},
                        {'Amino_Acid_Change': 'n.472G>T',
                         'Effect': 'non_coding_exon_variant',
                         'Effect_Impact': 'MODIFIER',
                         'Exon_Rank': 1,
                         'Gene_Coding': 'CODING',
                         'Gene_Name': 'RNF130',
                         'Genotype_Number': 1,
                         'Transcript_BioType': 'retained_intron',
                         'Transcript_ID': 'ENST00000521901'}],
                'MLEAC': '2',
                'MLEAF': '1.0'},
          'ALL': {'AN': '2',
                  'DP': '2',
                  'FS': '0.0',
                  'MQ': '60.0',
                  'MQ0': '0',
                  'QD': '10.88'}},
 'POS': 179390472,
 'QUAL': 21.77,
 'REF': 'C',
 'SAMPLES': {'TEST_CHR5': {'AD': [0, 2],
                           'DP': 2,
                           'GQ': 6,
                           'GT': '1/1',
                           'PL': [49, 6, 0]}}}
```
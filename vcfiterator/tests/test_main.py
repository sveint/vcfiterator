import unittest
from StringIO import StringIO

from vcfiterator import VcfIterator

class StringIOWrapper(StringIO):
    """
    Wrapper to support xreadlines()
    """
    def xreadlines(self):
        """
        We don't mind memory usage for these tests..
        """
        lines = self.readlines()
        for line in lines:
            yield line

HEADER = """##fileformat=VCFv4.1
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TESTSAMPLE1	TESTSAMPLE2	TESTSAMPLE3"""


def get_vcf_file_obj(variants):
    str_data = [HEADER]
    if variants:
        str_data.append(variants)
    return StringIOWrapper('\n'.join(str_data))


class TestHeaderParser(unittest.TestCase):

    def test_get_samples(self):
        vi = VcfIterator(get_vcf_file_obj(None))
        self.assertEquals(
            vi.getSamples(),
            ['TESTSAMPLE1', 'TESTSAMPLE2', 'TESTSAMPLE3']
        )

    def test_meta_parsing(self):
        vi = VcfIterator(get_vcf_file_obj(None))
        meta = vi.getMeta()
        self.assertEquals(
            meta['FORMAT'],
            [
                {
                    'Description': 'Genotype',
                     'Type': 'String',
                     'ID': 'GT',
                     'Number': '1'
                 },
                 {
                     'Description': 'Genotype Quality',
                     'Type': 'Integer',
                     'ID': 'GQ',
                     'Number': '1'
                 },
                 {
                     'Description': 'Read Depth',
                     'Type': 'Integer',
                     'ID': 'DP',
                     'Number': '1'
                 },
                 {
                    'Description': 'Haplotype Quality',
                    'Type': 'Integer',
                    'ID': 'HQ',
                    'Number': '2'
                }
            ]
        )


class TestDataParser(unittest.TestCase):

    def test_general_parsing(self):
        v = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.'
        vi = VcfIterator(get_vcf_file_obj(v))
        data = list(vi.iter())[0]

        self.assertEquals(data['CHROM'], '20')
        self.assertEquals(data['POS'], 14370)
        self.assertEquals(data['ID'], 'rs6054257')
        self.assertEquals(data['REF'], 'G')
        self.assertEquals(data['ALT'], ['A'])
        self.assertEquals(data['QUAL'], 29)
        self.assertEquals(data['FILTER'], 'PASS')

    def test_sample_parsing(self):
        v = '20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.'
        vi = VcfIterator(get_vcf_file_obj(v))
        data = list(vi.iter())[0]
        self.assertIn('TESTSAMPLE1', data['SAMPLES'])
        self.assertIn('TESTSAMPLE2', data['SAMPLES'])
        self.assertIn('TESTSAMPLE3', data['SAMPLES'])

        self.assertEquals(
            data['SAMPLES']['TESTSAMPLE1'],
            {
                'GT': '0|0',
                'GQ': 48,
                'DP': 1,
                'HQ': [51, 51]
            }
        )

        self.assertEquals(
            data['SAMPLES']['TESTSAMPLE2'],
            {
                'GT': '1|0',
                'GQ': 48,
                'DP': 8,
                'HQ': [51, 51]
            }
        )

        self.assertEquals(
            data['SAMPLES']['TESTSAMPLE3'],
            {
                'GT': '1/1',
                'GQ': 43,
                'DP': 5,
                'HQ': ['.', '.']
            }
        )

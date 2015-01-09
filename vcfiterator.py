import sys
from collections import defaultdict
import re


# Official fields in specification
SPEC_FIELDS = [
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER',
    'INFO',
    'FORMAT'
]


class Util(object):

    @staticmethod
    def conv_to_number(value):
        """
        Tries to convert a string to a number, silently returning the originally value if it fails.
        """
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        return value

    @staticmethod
    def split_and_convert(conv_func, split_max=-1, extract_single=False):
        """
        Performs a normal split() on a string, with support for converting the values and extraction of single values.

        :param conv_func: Function for converting the values
        :type conv_func: functions
        :param split_max: Maximum number of splits to perform. Default: No limit.
        :type split_max: int
        :param extract_single: If value ends up being a single value, do not return a list. Default: False
        :type extract_single: bool
        """
        def inner(x):
            l = [conv_func(i) for i in x.split(',', split_max)]
            if len(l) == 1 and extract_single:
                l = l[0]
            return l
        return inner


class VEPInfo(object):
    """
    Parser for the VEP INFO field.
    """

    field = 'CSQ'

    def __init__(self, meta):
        self.meta = meta
        self.fields = self._parseFieldsFromMeta()
        self.converters = {
            'DISTANCE': int,
            'ALLELE_NUM': int,
            'STRAND': int
        }

    def _parseFieldsFromMeta(self):
        info_line = next((l for l in self.meta['INFO'] if l.get('ID') == VEPInfo.field), None)
        if info_line:
            fields = info_line['Description'].split('Format: ', 1)[1].split('|')
            return fields
        return list()

    def parse(self, key, value, info_data, alleles):
        transcripts = value.split(',')

        all_data = [
            {
                k: self.converters.get(k, lambda x: x)(v) for k, v in zip(self.fields, t.split('|')) if v is not ''
            } for t in transcripts
        ]

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = [d for d in all_data if d['ALLELE_NUM']-1 == a_idx]


class SnpEffInfo(object):
    """
    Parser for the snpEff INFO field.
    """

    field = 'EFF'

    def __init__(self, meta):
        self.meta = meta
        self.fields = self._parseFieldsFromMeta()
        self.converters = {
            'Genotype_Number': int,
            'Exon_Rank': int,
            'Amino_Acid_length': int
        }

    def _parseFormat(self, line):
        """
        Parse following format to a flat list.
        Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )
        """
        fields = list()

        line = line.replace('(', '|').replace(')', '').replace('[ | ERRORS | WARNINGS ]', '').replace('\'', '')
        fields = line.split('|')

        fields = [f.strip() for f in fields]

        return fields

    def _parseFieldsFromMeta(self):
        info_line = next((l for l in self.meta['INFO'] if l.get('ID') == SnpEffInfo.field), None)
        if info_line:
            fields = self._parseFormat(info_line['Description'].split('Format: \'', 1)[1])
            fields.append('ERRORS')
            return fields
        return list()

    def parse(self, key, value, info_data, alleles):
        transcripts = value.split(',')

        all_data = [
            {
                k: self.converters.get(k, lambda x: x)(v) for k, v in zip(self.fields, self._parseFormat(t)) if v is not ''
            } for t in transcripts
        ]

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = [d for d in all_data if d['Genotype_Number']-1 == a_idx]


class CsvAlleleParser(object):
    """
    Parses comma separated values, and inserts them into the data according to the allele the value belongs to.
    """
    def __init__(self, meta):
        self.meta = meta

    def parse(self, key, value, info_data, alleles):
        allele_values = value.split(',')
        if not len(allele_values) == len(alleles):
            raise RuntimeError("Number of allele values for {} not matching number of alleles".format(key))

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = allele_values[a_idx]


class HeaderParser(object):
    """
    Class for parsing the header part of the vcf and returns the metadata and header data.
    """

    RE_INFO = re.compile(r'[<]*(.*?)=["]*(.*?)["]*[,>]')

    def __init__(self, path):
        self.path = path
        self.metaProccessors = {
            'INFO': self._parseMetaInfo,
            'FILTER': self._parseMetaInfo,
            'FORMAT': self._parseMetaInfo
        }

    def _getSamples(self, header):
        return [field for field in header if field not in SPEC_FIELDS]

    def _parseMetaInfo(self, infoline):
        groups = re.findall(HeaderParser.RE_INFO, infoline)
        info = {k: v for k, v in groups}
        return info

    def _parseHeader(self):
        meta = defaultdict(list)
        header = list()

        # Read in metadata and header
        with open(self.path) as fd:
            for line in fd.xreadlines():
                line = line.replace('\n', '')
                if line.startswith('##'):
                    key, value = line[2:].split('=', 1)
                    meta[key].append(value)
                elif(line.startswith('#')):
                    line = line.replace('#', '')
                    header = line.split('\t')
                else:
                    # End of header
                    break

        # Extract data with processors
        for key, func in self.metaProccessors.iteritems():
            if key in meta:
                for idx, value in enumerate(meta[key]):
                    meta[key][idx] = func(value)

        # Extract value from single-item lists ([val] -> val):
        for k, v in meta.iteritems():
            if len(v) == 1:
                meta[k] = v[0]

        samples = self._getSamples(header)
        return meta, header, samples

    def parse(self):
        return self._parseHeader()


class DataParser(object):

    def __init__(self, path, meta, header, samples):
        self.path = path
        self.meta = meta
        self.header = header
        self.samples = samples

        self.infoProcessors = self._generateInfoProcessors()
        self.infoProcessors.update({
            SnpEffInfo.field: SnpEffInfo(self.meta).parse,
            VEPInfo.field: VEPInfo(self.meta).parse,
            'AC': CsvAlleleParser(self.meta).parse,
            'AF': CsvAlleleParser(self.meta).parse,
            'MLEAC': CsvAlleleParser(self.meta).parse,
            'MLEAF': CsvAlleleParser(self.meta).parse,
        })

    def _generateInfoProcessors(self):
        """
        Generates the different processors to handle the keys in the INFO field.

        These processors are generated from INFO fields in the header metadata, using the specified type and length.

        Unlike the custom processors (e.g. VEPInfo), it does not generate allele specific data. Instead, all data is put
        into the 'ALL' key in 'INFO' in the resulting dictionary.
        """
        processors = dict()
        for f in self.meta['INFO']:
            parse_func = str
            if f['Type'] == 'Integer':
                parse_func = int
            elif f['Type'] in ['Number', 'Double', 'Float']:
                parse_func = float
            elif f['Type'] == 'Flag':
                parse_func = bool

            number = f['Number']

            try:
                # Number == int
                n = int(number)
                func = Util.split_and_convert(parse_func, split_max=n, extract_single=True)
            except ValueError:
                # Number == Allele specific
                if number == 'A':
                    func = Util.split_and_convert(parse_func)
                # Number == Unknown
                else:
                    func = parse_func

            # We ignore alleles for these values, but return them in the 'ALL' key
            def val_func(key, value, info_data, alleles):
                info_data['ALL'][key] = func(value)

            # Add processor for this value, recognized by it's ID
            processors[f['ID']] = val_func

        return processors

    def _parseDataInfoField(self, data):
        """
        Parses the INFO data into data structures.
        Data is split into general ('ALL') and allele specific data.
        """

        alleles = data['ALT']

        fields = data['INFO'].split(';')

        # Create dict for allele specific INFO
        info_data = {
            k: dict() for k in alleles
        }
        # And include INFO for 'ALL' alleles
        info_data['ALL'] = dict()

        for f in fields:
            if '=' in f:
                key, value = f.split('=', 1)
                # Process keys by processor, if present, or just give value
                # Data is inserted into info_data by the functions
                self.infoProcessors[key](key, value, info_data, alleles)
            else:
                info_data['ALL'][f] = True

        data['INFO'] = info_data

    def _parseDataSampleFields(self, data):
        sample_format = data['FORMAT'].split(':')

        samples = dict()
        extract = Util.split_and_convert(Util.conv_to_number, extract_single=True)
        for sample_name in self.samples:
            sample_text = data.pop(sample_name)
            samples[sample_name] = {
                k: extract(v) for k, v in zip(sample_format, sample_text.split(':'))
            }

        data['SAMPLES'] = samples

        del data['FORMAT']

    def _parseData(self, line):
        data = {
            k: v for k, v in zip(self.header, line.split('\t'))
        }

        # Split by alleles
        data['ALT'] = data['ALT'].split(',')

        self._parseDataInfoField(data)

        self._parseDataSampleFields(data)

        # Manual conversion
        data['POS'] = int(data['POS'])
        data['QUAL'] = float(data['QUAL'])

        return data

    def iter(self):
        found_data_start = False
        with open(self.path) as fd:
            for line in fd.xreadlines():
                # Skip header, wait for #CHROM to signal start of data
                if line.startswith('#CHROM') and not found_data_start:
                    found_data_start = True
                    continue
                if not found_data_start:
                    continue
                line = line.replace('\n', '')
                data = self._parseData(line)
                yield data


class VcfIterator(object):

    def __init__(self, path):
        self.path = path

        self.meta, self.header, self.samples = HeaderParser(self.path).parse()

    def getHeader(self):
        return self.header

    def getMeta(self):
        return self.meta

    def getSamples(self):
        return self.samples

    def iter(self):
        d = DataParser(self.path, self.meta, self.header, self.samples)
        for r in d.iter():
            yield r


if __name__ == '__main__':
    import json

    path = sys.argv[1]
    v = VcfIterator(path)

    for value in v.iter():
        print json.dumps(value, indent=4)


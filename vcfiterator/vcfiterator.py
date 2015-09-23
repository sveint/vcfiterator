import sys
from collections import defaultdict
import re
import abc


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


class BaseInfoProcessor(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def accepts(self, key, value, processed):
        """
        Checks whether this info processor should be run for this key/value.

        :param key: The key of the INFO field.
        :param value: The string value for this field.
        :param processed: Tells whether another processor has already accepted this field.

        """
        pass

    @abc.abstractmethod
    def process(self, key, value, info_data, alleles, processed):
        """
        For processing the incoming key, value pair, inserting the data into info_data
        however is seen fit. Is only invoked if accepts returned True.

        :param key: The key of the INFO field.
        :param value: The string value for this field.
        :param info_data: INFO data structure for inserting data into.
        :param alleles: List of alleles (strings) for this value. In practice same as ALT field.
        :param processed: Tells whether another processor has already accepted this field.

        """
        pass

    def getConvertFunction(self, meta, key):
        # Search for meta item
        f = next((m for m in meta['INFO'] if m['ID'] == key), None)
        func = lambda x: x.decode('latin-1', 'replace')
        if f:
            parse_func = str
            if f['Type'] == 'Integer':
                parse_func = int
            elif f['Type'] in ['Number', 'Double', 'Float']:
                parse_func = float
            elif f['Type'] == 'Flag':
                parse_func = bool
            elif f['Type'] == 'String':
                parse_func = lambda x: x.decode('latin-1', 'replace')

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
                    if f['Type'] == 'Integer':
                        func = Util.split_and_convert(parse_func, extract_single=True)
                    else:
                        func = parse_func

        return func


class VEPInfoProcessor(BaseInfoProcessor):
    """
    Parser for the VEP INFO field.
    """

    field = 'CSQ'

    def __init__(self, meta):
        self.meta = meta
        self.fields = self._parseFieldsFromMeta()
        self.converters = {
            'AA_MAF': self._parseMAF,
            'AFR_MAF': self._parseMAF,
            'AMR_MAF': self._parseMAF,
            'ALLELE_NUM': int,
            'ASN_MAF': self._parseMAF,
            'EA_MAF': self._parseMAF,
            'EUR_MAF': self._parseMAF,
            'EAS_MAF': self._parseMAF,
            'SAS_MAF': self._parseMAF,
            'GMAF': self._parseMAF,
            'EAS_MAF': self._parseMAF,
            'SAS_MAF': self._parseMAF,
            'Consequence': lambda x: [i for i in x.split('&')],
            'Existing_variation': lambda x: [i for i in x.split('&')],
            'DISTANCE': int,
            'STRAND': int,
            'PUBMED': lambda x: [int(i) for i in x.split('&')],
        }

    def _parseFieldsFromMeta(self):
        info_line = next((l for l in self.meta['INFO'] if l.get('ID') == VEPInfoProcessor.field), None)
        if info_line:
            fields = info_line['Description'].split('Format: ', 1)[1].split('|')
            return fields
        return list()

    def _parseMAF(self, val):
        maf = dict()
        alleles = val.split('&')
        for allele in alleles:
            v = allele.split(':')
            for key, value in zip(v[0::2], v[1::2]):
                try:
                    maf[key] = float(value)
                except ValueError:
                    continue
        return maf

    def accepts(self, key, value, processed):
        return key == VEPInfoProcessor.field

    def process(self, key, value, info_data, alleles, processed):
        transcripts = value.split(',')

        all_data = [
            {
                k: self.converters.get(k, lambda x: x.decode('latin-1', 'replace'))(v) for k, v in zip(self.fields, t.split('|')) if v is not ''
            } for t in transcripts
        ]

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = [d for d in all_data if d['ALLELE_NUM']-1 == a_idx]


class SnpEffInfoProcessor(BaseInfoProcessor):
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
        info_line = next((l for l in self.meta['INFO'] if l.get('ID') == SnpEffInfoProcessor.field), None)
        if info_line:
            fields = self._parseFormat(info_line['Description'].split('Format: \'', 1)[1])
            fields.append('ERRORS')
            return fields
        return list()

    def accepts(self, key, value, processed):
        return key == SnpEffInfoProcessor.field

    def process(self, key, value, info_data, alleles, processed):
        transcripts = value.split(',')

        all_data = [
            {
                k: self.converters.get(k, lambda x: x.decode('latin-1', 'replace'))(v) for k, v in zip(self.fields, self._parseFormat(t)) if v is not ''
            } for t in transcripts
        ]

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = [d for d in all_data if d['Genotype_Number']-1 == a_idx]


class CsvAlleleParser(BaseInfoProcessor):
    """
    Parses comma separated values, and inserts them into the data according to the allele the value belongs to.
    """

    fields = ['AC', 'AF', 'MLEAC', 'MLEAF']

    def __init__(self, meta):
        self.meta = meta
        self.conv_func = Util.conv_to_number

    def accepts(self, key, value, processed):
        return key in CsvAlleleParser.fields

    def process(self, key, value, info_data, alleles, processed):
        allele_values = value.split(',')
        if not len(allele_values) == len(alleles):
            raise RuntimeError("Number of allele values for {} not matching number of alleles".format(key))

        for a_idx, allele in enumerate(alleles):
            info_data[allele][key] = self.conv_func(allele_values[a_idx])


class NativeInfoProcessor(BaseInfoProcessor):
        """
        Fallback processor, invoked if none of the custom ones accepted the data.

        It searches the INFO fields in the header metadata, trying to use the specified type and length.

        Unlike custom processors (e.g. VEPInfoProcessor), it does not generate allele specific data. Instead, all data is put
        into the 'ALL' key in 'INFO' in the resulting dictionary.
        """

        def __init__(self, meta):
            self.meta = meta

        def accepts(self, key, value, processed):
            return not processed

        def process(self, key, value, info_data, alleles):

            if isinstance(value, bool):
                info_data['ALL'][key] = value
            else:
                func = self.getConvertFunction(self.meta, key)
                # We ignore alleles for these values, but return them in the 'ALL' key
                info_data['ALL'][key] = func(value)


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

        self.infoProcessors = list()
        self.fallbackProcessor = NativeInfoProcessor(meta)

    def addInfoProcessor(self, processor):
        self.infoProcessors.append(processor)

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
            else:
                key, value = f, True
            # Process keys by processor, if present, or use native processor
            # Data is inserted into info_data by the functions
            processed = False
            for processor in self.infoProcessors:
                if processor.accepts(key, value, processed):
                    processor.process(key, value, info_data, alleles, processed)
                    processed = True
            # If no processors handled the data, use the native header processor
            if not processed:
                self.fallbackProcessor.process(key, value, info_data, alleles)

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
        data['POS'] = Util.conv_to_number(data['POS'])
        data['QUAL'] = Util.conv_to_number(data['QUAL'])

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
        self.data_parser = DataParser(self.path, self.meta, self.header, self.samples)

        self.addInfoProcessor(VEPInfoProcessor(self.meta))
        self.addInfoProcessor(SnpEffInfoProcessor(self.meta))
        self.addInfoProcessor(CsvAlleleParser(self.meta))

    def getHeader(self):
        return self.header

    def getMeta(self):
        return self.meta

    def getSamples(self):
        return self.samples

    def addInfoProcessor(self, processor):
        self.data_parser.addInfoProcessor(processor)

    def iter(self):
        for r in self.data_parser.iter():
            yield r


if __name__ == '__main__':
    import json

    path = sys.argv[1]
    v = VcfIterator(path)

    for value in v.iter():
        print json.dumps(value, indent=4)


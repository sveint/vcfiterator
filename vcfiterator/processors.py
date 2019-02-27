import abc

from vcfiterator.util import Util


class BaseInfoProcessor(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, meta):
        self.meta = meta

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
            parse_func = Util.dot_to_none(lambda x: x.decode('latin-1', 'replace'))
            if f['Type'] == 'Integer':
                parse_func = Util.dot_to_none(int)
            elif f['Type'] in ['Number', 'Double', 'Float']:
                parse_func = Util.dot_to_none(float)
            elif f['Type'] == 'Flag':
                parse_func = Util.dot_to_none(bool)

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
        super(VEPInfoProcessor, self).__init__(meta)

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

        if len(alleles) == 1:
            info_data[alleles[0]][key] = all_data
        else:
            for a_idx, allele in enumerate(alleles):
                info_data[allele][key] = [d for d in all_data if d['ALLELE_NUM']-1 == a_idx]


class SnpEffInfoProcessor(BaseInfoProcessor):
    """
    Parser for the snpEff INFO field.
    """

    field = 'EFF'

    def __init__(self, meta):
        super(SnpEffInfoProcessor, self).__init__(meta)

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
        super(CsvAlleleParser, self).__init__(meta)
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
            super(NativeInfoProcessor, self).__init__(meta)

        def accepts(self, key, value, processed):
            return not processed

        def process(self, key, value, info_data, alleles):

            if isinstance(value, bool):
                info_data['ALL'][key] = value
            else:
                func = self.getConvertFunction(self.meta, key)
                # We ignore alleles for these values, but return them in the 'ALL' key
                info_data['ALL'][key] = func(value)

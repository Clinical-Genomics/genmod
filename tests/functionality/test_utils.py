from typing import Dict, Union
from genmod.vcf_tools import HeaderParser, get_variant_dict, get_info_dict

def parse_variant_file(file_path: str) -> HeaderParser:
    """
    Parse VCF header fields
    :param file_path: VCF to be read
    :raises ValueError: in case file is empty
    """
    with open(file_path, 'r') as variant_file:
        head = HeaderParser()
        for line_index, line in enumerate(variant_file):
            line = line.rstrip()
            if line.startswith('#'):
                if line.startswith('##'):
                    head.parse_meta_data(line)
                else:
                    head.parse_header_line(line)
            else:
                break
        if line_index == 0:
            raise ValueError('Expected contents in file, got none')
    return head

def generate_variants_from_file(file_path: str) -> Dict[str, Union[str, int, float]]:
    """
    Yield variants from VCF file.
    :param file_path: VCF to be read
    """
    header = parse_variant_file(file_path=file_path)
    with open(file_path, 'r') as variant_file:
        for line in variant_file:
            if line.startswith('#'):
                continue
            variant: Dict[str, str] = get_variant_dict(line, header.header)
            variant['info_dict'] = get_info_dict(variant['INFO'])
            yield variant


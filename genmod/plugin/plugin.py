"""Vendored from extract_vcf (https://github.com/moonso/extract_vcf, MIT License, Måns Magnusson)."""

from __future__ import annotations

import operator
from collections.abc import Sequence
from logging import getLogger
from typing import Any, Optional, Union

from .split_strings import split_strings

VariantDict = dict[str, Any]
"""A variant as a dictionary, with VCF columns as keys and ``info_dict`` for INFO."""

Value = Union[str, int, float, bool, None]
"""The value types returned by :meth:`Plugin.get_value`."""


class Plugin(object):
    """Extract a single value from a VCF variant based on a set of rules.

    A plugin describes where in a variant to look (``field``, ``info_key``, ``csq_key``,
    ``gt_key``), how the record is split up (``separators``), what type it has
    (``data_type``) and how to choose between multiple values (``record_rule`` and
    ``string_rules``). Use :meth:`get_value` to apply it to a variant, given either as a
    VCF line or as a dictionary.

    Attributes:
        data_types: The valid values for ``data_type``.
        regex: Unused, kept for compatibility with ``extract_vcf``.
    """

    def __init__(
        self,
        name: str,
        field: str,
        data_type: Optional[str] = None,
        separators: Sequence[str] = [],
        info_key: Optional[str] = None,
        category: Optional[str] = None,
        csq_key: Optional[str] = None,
        record_rule: Optional[str] = None,
        gt_key: Optional[str] = None,
        string_rules: dict[str, int] = {},
        dict_entry: bool = False,
    ) -> None:
        """Initialise the plugin.

        Args:
            name: The name of the plugin.
            field: The field of the VCF that should be searched. Any of
                ``CHROM``, ``POS``, ``ID``, ``REF``, ``ALT``, ``QUAL``, ``FILTER``,
                ``INFO``, ``FORMAT`` or ``sample_id``.
            data_type: The data type of the record. Any of ``integer``, ``float``,
                ``flag`` or ``string``.
            separators: How the record is separated, applied in order.
            info_key: The name of the INFO field (when ``field`` is ``INFO``).
            category: The category the plugin belongs to.
            csq_key: The name of the VEP entry (when ``info_key`` is ``CSQ``).
            record_rule: Which value to return if there are several, ``min`` or ``max``.
            gt_key: The FORMAT key to extract (when ``field`` is ``sample_id``).
            string_rules: Priority of string matches, mapping the string to its rank.
                Used with ``data_type="string"`` and a ``record_rule``.
            dict_entry: If the values are annotated as a dictionary. The first two
                separators then split the dictionary and its keys and values, and the
                remaining ones split the value.
        """
        super(Plugin, self).__init__()
        self.logger = getLogger(__name__)
        # For testing only
        self.logger = getLogger("genmod.plugin")

        # These are the valid data types for a entry
        self.data_types = ["integer", "float", "flag", "string"]
        self.name = name
        self.logger.info("Initiating plugin with name: {0}".format(self.name))
        self.regex = None
        self.field = field
        self.logger.info("Field: {0}".format(self.field))

        self.data_type = data_type
        self.logger.info("Data type: {0}".format(self.data_type))

        self.separators = separators
        self.logger.info("Separators: {0}".format(self.separators))

        self.record_rule = record_rule
        self.logger.info("Record rule: {0}".format(self.record_rule))

        self.info_key = info_key
        self.logger.info("Info key: {0}".format(self.info_key))

        self.csq_key = csq_key
        self.logger.info("CSQ key: {0}".format(self.csq_key))

        self.category = category
        self.logger.info("Category: {0}".format(self.category))

        self.string_rules = string_rules
        self.logger.info("String rules: {0}".format(self.string_rules))

        self.gt_key = gt_key
        self.logger.info("gt_key: {0}".format(self.gt_key))

        self.dict_entry = dict_entry

    def get_entry(
        self,
        variant_line: Optional[str] = None,
        variant_dict: Optional[VariantDict] = None,
        raw_entry: Optional[str] = None,
        vcf_header: Optional[Sequence[str]] = None,
        csq_format: Optional[Sequence[str]] = None,
        dict_key: Optional[str] = None,
        individual_id: Optional[str] = None,
    ) -> list[str]:
        """Return the entry of the variant, split into its separate values.

        Args:
            variant_line: A VCF formatted variant line.
            variant_dict: A variant dictionary, used if no ``variant_line`` is given.
            raw_entry: The raw entry. If given it is split directly instead of being
                looked up in the variant.
            vcf_header: The VCF header line (column names including sample ids).
            csq_format: The names of the CSQ (VEP) columns. Required for ``CSQ`` entries.
            dict_key: The key to pick if the entry is a dictionary (see ``dict_entry``).
            individual_id: The sample id to use when ``field`` is ``sample_id``.

        Returns:
            The split entry. Empty if no entry was found.

        Raises:
            IOError: If the CSQ format or CSQ key is missing for a CSQ entry.
        """
        if not raw_entry:
            raw_entry = self.get_raw_entry(
                variant_line=variant_line,
                variant_dict=variant_dict,
                vcf_header=vcf_header,
                individual_id=individual_id,
                dict_key=dict_key,
            )

        entry = []
        if raw_entry:
            if self.field in ["CHROM", "POS", "REF", "QUAL"]:
                # We know these fields allways has one entry
                entry = [raw_entry]

            elif self.field in ["ID", "FILTER"]:
                # We know ID is allways splitted on ';'
                entry = raw_entry.split(";")

            elif self.field == "ALT":
                # We know ALT is allways splitted on ','
                entry = raw_entry.split(",")

            elif self.field == "FORMAT":
                entry = raw_entry.split(":")

            elif self.field == "INFO":
                # We are going to treat csq fields separately
                if self.info_key == "CSQ":
                    if not csq_format:
                        raise IOError("If CSQ the csq format must be provided")
                    if not self.csq_key:
                        raise IOError("If CSQ a csq key must be provided")
                    for i, head in enumerate(csq_format):
                        if head == self.csq_key:
                            # This is the csq entry we are looking for
                            csq_column = i
                    # CSQ entries are allways splitted on ','
                    for csq_entry in raw_entry.split(","):
                        entry += split_strings(csq_entry.split("|")[csq_column], self.separators)
                else:
                    if self.dict_entry:
                        separators = self.separators[2:]
                    else:
                        separators = self.separators

                    entry = split_strings(raw_entry, separators)

            elif self.field == "sample_id":
                if not self.separators:
                    entry = split_strings(raw_entry, "/")
                    # If variant calls are phased we need to split on '|'
                    if len(entry) == 1:
                        entry = split_strings(raw_entry, "|")
                else:
                    entry = split_strings(raw_entry, self.separators)

        return entry

    def get_raw_entry(
        self,
        variant_line: Optional[str] = None,
        variant_dict: Optional[VariantDict] = None,
        vcf_header: Optional[Sequence[str]] = None,
        individual_id: Optional[str] = None,
        dict_key: Optional[str] = None,
    ) -> Optional[str]:
        """Return the raw, unsplit entry from the VCF field of the plugin.

        Args:
            variant_line: A VCF formatted variant line.
            variant_dict: A variant dictionary, used if no ``variant_line`` is given.
            vcf_header: The VCF header line. Required for ``sample_id`` with a line.
            individual_id: The sample id to get the genotype call for.
            dict_key: The key to pick if the entry is a dictionary (see ``dict_entry``).

        Returns:
            The raw entry, or None if no entry was found.

        Raises:
            IOError: If ``field`` is ``sample_id`` and the individual id, genotype key
                or (for a variant line) the VCF header is missing.
        """
        columns = variant_line.rstrip().split() if variant_line else None

        entry = None

        if self.field == "CHROM":
            if columns:
                entry = columns[0]
            elif variant_dict:
                entry = variant_dict["CHROM"]

        elif self.field == "POS":
            if columns:
                entry = columns[1]
            elif variant_dict:
                entry = variant_dict["POS"]

        elif self.field == "ID":
            if columns:
                entry = columns[2]
            elif variant_dict:
                entry = variant_dict["ID"]

        elif self.field == "REF":
            if columns:
                entry = columns[3]
            elif variant_dict:
                entry = variant_dict["REF"]

        elif self.field == "ALT":
            if columns:
                entry = columns[4]
            elif variant_dict:
                entry = variant_dict["ALT"]

        elif self.field == "QUAL":
            if columns:
                entry = columns[5]
            elif variant_dict:
                entry = variant_dict["QUAL"]

        elif self.field == "FILTER":
            if columns:
                entry = columns[6]
            elif variant_dict:
                entry = variant_dict["FILTER"]

        elif self.field == "INFO":
            if columns:
                for info_annotation in columns[7].split(";"):
                    splitted_annotation = info_annotation.split("=")
                    if self.info_key == splitted_annotation[0]:
                        if len(splitted_annotation) == 2:
                            entry = splitted_annotation[1]

            elif variant_dict:
                entry = variant_dict.get("info_dict", {}).get(self.info_key)

            if self.dict_entry and entry:
                # First we split the "dictionaries"
                first_split = entry.split(self.separators[0])
                for annotation in first_split:
                    # Then we search for the dict key
                    splitted_entry = annotation.split(self.separators[1])
                    key = splitted_entry[0]
                    value = splitted_entry[1]
                    if dict_key:
                        if key == dict_key:
                            entry = value
                    # If no key we just return the last entry
                    else:
                        entry = value

        elif self.field == "FORMAT":
            if columns:
                entry = columns[8]
            elif variant_dict:
                entry = variant_dict["FORMAT"]

        elif self.field == "sample_id":
            if not individual_id:
                raise IOError("If 'sample_id' a individual id must be provided")
            if not self.gt_key:
                raise IOError("If 'sample_id' a genotype key must be provided")

            if columns:
                if not vcf_header:
                    raise IOError("If 'sample_id' the vcf header must be provided")

                format_info = columns[8]

                for i, head in enumerate(vcf_header):
                    if head == individual_id:
                        raw_gt_call = columns[i]
            elif variant_dict:
                format_info = variant_dict["FORMAT"]
                raw_gt_call = variant_dict[individual_id]

            entry_dict = dict(zip(format_info.split(":"), raw_gt_call.split(":")))
            entry = entry_dict.get(self.gt_key, ".")

        return entry

    def get_value(
        self,
        variant_line: Optional[str] = None,
        variant_dict: Optional[VariantDict] = None,
        entry: Optional[Sequence[str]] = None,
        raw_entry: Optional[str] = None,
        vcf_header: Optional[Sequence[str]] = None,
        csq_format: Optional[Sequence[str]] = None,
        dict_key: Optional[str] = None,
        individual_id: Optional[str] = None,
    ) -> Value:
        """Return the value of the variant as specified by the plugin.

        ``entry`` and ``raw_entry`` are accepted for compatibility with ``extract_vcf``
        but are not used; the raw entry is always looked up from the variant.

        Args:
            variant_line: A VCF formatted variant line.
            variant_dict: A variant dictionary, used if no ``variant_line`` is given.
            entry: Unused.
            raw_entry: Unused.
            vcf_header: The VCF header line with sample ids.
            csq_format: The names of the CSQ (VEP) columns.
            dict_key: The key to pick if the entry is a dictionary (see ``dict_entry``).
            individual_id: The sample id to use when ``field`` is ``sample_id``.

        Returns:
            The value, typed according to ``data_type``: ``True`` for a ``flag`` that is
            present, the matching string rule for a ``string``, or the ``int``/``float``
            chosen by ``record_rule``. None if no value was found.
        """
        # Also used for intermediate values while looping, so not just ``Value``
        value: Any = None

        raw_entry = self.get_raw_entry(
            variant_line=variant_line,
            variant_dict=variant_dict,
            vcf_header=vcf_header,
            individual_id=individual_id,
            dict_key=dict_key,
        )
        # If data type is flag we only need to check if any entry exists
        if self.data_type == "flag":
            if self.field == "INFO":
                if variant_line:
                    for info_entry in variant_line.split()[7].split(";"):
                        if self.info_key == info_entry.split("=")[0]:
                            value = True
                elif variant_dict:
                    if self.info_key in variant_dict.get("info_dict", {}):
                        value = True
            else:
                if raw_entry != ".":
                    value = True

        # If we have a record rule we need to return the correct value
        elif raw_entry:
            # If there was no raw entry we will return None
            if self.record_rule:
                if self.data_type == "string":
                    if self.record_rule == "max":
                        sorted_strings = sorted(
                            self.string_rules.items(), key=operator.itemgetter(1), reverse=True
                        )

                    if self.record_rule == "min":
                        sorted_strings = sorted(
                            self.string_rules.items(), key=operator.itemgetter(1)
                        )

                    for string_rule in sorted_strings:
                        if string_rule[0].lower() in raw_entry.lower():
                            value = string_rule[0]
                            break
                else:
                    typed_annotations = []

                    for value in self.get_entry(
                        raw_entry=raw_entry,
                        vcf_header=vcf_header,
                        csq_format=csq_format,
                        dict_key=dict_key,
                        individual_id=individual_id,
                    ):
                        if self.data_type == "float":
                            try:
                                typed_annotations.append(float(value))
                            except ValueError:
                                pass

                        elif self.data_type == "integer":
                            try:
                                typed_annotations.append(int(value))
                            except ValueError:
                                pass

                    if typed_annotations:
                        if self.record_rule == "max":
                            value = max(typed_annotations)

                        elif self.record_rule == "min":
                            value = min(typed_annotations)
                    else:
                        value = None

            # If no record rule is given we return the raw annotation
            # Here the data_type is not flag, and there is no record rule
            # We know that there exists a raw annotation
            else:
                # We will just return the first annotation found
                value = self.get_entry(
                    raw_entry=raw_entry,
                    vcf_header=vcf_header,
                    csq_format=csq_format,
                    dict_key=dict_key,
                    individual_id=individual_id,
                )[0]

                if self.data_type == "float":
                    try:
                        value = float(value)
                    except ValueError:
                        pass

                elif self.data_type == "integer":
                    try:
                        value = int(value)
                    except ValueError:
                        pass

        return value

    def __repr__(self) -> str:
        """Return a string with the main settings of the plugin."""
        return (
            "Plugin(name={0},field={1},data_type={2},separators={3},"
            "record_rule={4},info_key={5},csq_key={6},category={7},"
            "string_rules={8})".format(
                self.name,
                self.field,
                self.data_type,
                self.separators,
                self.record_rule,
                self.info_key,
                self.csq_key,
                self.category,
                self.string_rules,
            )
        )

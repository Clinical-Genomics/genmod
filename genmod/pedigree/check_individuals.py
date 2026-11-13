"""Check a pedigree against the samples of a VCF file."""

from __future__ import annotations

from collections.abc import Collection, Iterable


def check_individuals_in_vcf(
    ped_individuals: Iterable[str], vcf_individuals: Collection[str]
) -> None:
    """Check that all individuals in the ped file are samples in the VCF file.

    Individuals that are in the VCF but not in the ped file are allowed.

    Args:
        ped_individuals: The individual ids of the ped file, e.g. ``FamilyParser.individuals``
        vcf_individuals: The sample names of the VCF file

    Raises:
        IOError: If any ped individual is missing from the VCF, naming the missing ones
    """
    missing = [individual for individual in ped_individuals if individual not in vcf_individuals]
    if missing:
        raise IOError(
            f"Individuals in PED file must exist in VCF file, missing: {', '.join(missing)}"
        )

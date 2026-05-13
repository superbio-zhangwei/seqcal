"""seqcal package."""

from .core import (
    extract_region_sequence,
    find_target_regions,
    read_fasta,
    reverse_complement,
    run_seqcal,
)

__all__ = [
    "extract_region_sequence",
    "find_target_regions",
    "read_fasta",
    "reverse_complement",
    "run_seqcal",
]

__version__ = "0.1.0"

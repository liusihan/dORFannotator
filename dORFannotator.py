# -*- coding: utf-8 -*-
"""
dORFannotator - High-performance 3'UTR variant annotation tool for downstream ORFs

@author: Sihan Liu
@version: 1.0.0
@date: 2025
@institution: Institute of Rare Diseases / West China Hospital
@license: GNU General Public License v3
"""

import os
import sys
import re
import gzip
import time
import logging
import argparse
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field
from collections import defaultdict, OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import multiprocessing as mp

# =============================================================================
# Dependency Checking
# =============================================================================

# Required packages with minimum versions
REQUIRED_PACKAGES = {
    'pysam': '0.19.0',
    'cyvcf2': '0.30.0',
    'biopython': '1.79',
    'pandas': '1.3.0',
    'numpy': '1.20.0',
}

def check_dependencies() -> bool:
    """
    Check all required dependencies and their versions.
    
    Returns:
        True if all dependencies are satisfied, False otherwise.
    """
    missing = []
    version_errors = []
    
    # Check pysam
    try:
        import pysam
        version = pysam.__version__
        if not _check_version(version, REQUIRED_PACKAGES['pysam']):
            version_errors.append(f"pysam: found {version}, required >= {REQUIRED_PACKAGES['pysam']}")
    except ImportError:
        missing.append('pysam')
    
    # Check cyvcf2
    try:
        import cyvcf2
        version = cyvcf2.__version__
        if not _check_version(version, REQUIRED_PACKAGES['cyvcf2']):
            version_errors.append(f"cyvcf2: found {version}, required >= {REQUIRED_PACKAGES['cyvcf2']}")
    except ImportError:
        missing.append('cyvcf2')
    
    # Check biopython
    try:
        import Bio
        version = Bio.__version__
        if not _check_version(version, REQUIRED_PACKAGES['biopython']):
            version_errors.append(f"biopython: found {version}, required >= {REQUIRED_PACKAGES['biopython']}")
    except ImportError:
        missing.append('biopython')
    
    # Check pandas
    try:
        import pandas
        version = pandas.__version__
        if not _check_version(version, REQUIRED_PACKAGES['pandas']):
            version_errors.append(f"pandas: found {version}, required >= {REQUIRED_PACKAGES['pandas']}")
    except ImportError:
        missing.append('pandas')
    
    # Check numpy
    try:
        import numpy
        version = numpy.__version__
        if not _check_version(version, REQUIRED_PACKAGES['numpy']):
            version_errors.append(f"numpy: found {version}, required >= {REQUIRED_PACKAGES['numpy']}")
    except ImportError:
        missing.append('numpy')
    
    # Report errors
    if missing:
        print(f"\nError: Missing required packages: {', '.join(missing)}")
        print(f"Please install with: pip install {' '.join(missing)}")
        return False
    
    if version_errors:
        print("\nWarning: Package version issues detected:")
        for err in version_errors:
            print(f"  - {err}")
        print("Consider upgrading packages for best compatibility.")
        # Continue with warning, don't fail
    
    return True


def _check_version(installed: str, required: str) -> bool:
    """Compare version strings (simple major.minor.patch comparison)."""
    try:
        installed_parts = [int(x) for x in installed.split('.')[:3]]
        required_parts = [int(x) for x in required.split('.')[:3]]
        
        # Pad with zeros
        while len(installed_parts) < 3:
            installed_parts.append(0)
        while len(required_parts) < 3:
            required_parts.append(0)
        
        return installed_parts >= required_parts
    except (ValueError, AttributeError):
        return True  # If can't parse, assume OK


# Run dependency check before importing
if not check_dependencies():
    sys.exit(1)

# Now import dependencies
import pysam
from cyvcf2 import VCF
from Bio.Seq import Seq, reverse_complement, translate
from Bio.SeqUtils import gc_fraction
import pandas as pd
import numpy as np

# =============================================================================
# Constants and Configuration
# =============================================================================

__version__ = '1.0.0'

# Default biological parameters (user configurable)
DEFAULT_START_CODONS = {'ATG', 'CTG', 'GTG', 'TTG', 'ACG'}
DEFAULT_STOP_CODONS = {'TAG', 'TAA', 'TGA'}
MIN_ORF_LENGTH = 30   # 10 aa
MAX_ORF_LENGTH = 303  # 101 aa

# Variant classification types
VARIANT_TYPES = [
    'dStart_gained',    # New start codon created
    'dStop_gained',     # New stop codon created
    'dStart_lost',      # Existing start codon destroyed
    'dStop_lost',       # Existing stop codon destroyed
    'dStart_change',    # Start codon changed to another valid start
    'dStop_change',     # Stop codon changed to another valid stop
    'dMissense',        # Missense variant in dORF
    'dIndel_frameshift', # Indel causing frameshift
    'dIndel_inframe'    # Indel maintaining reading frame
]

# Output column definitions
# Column name -> VariantEffect attribute mapping
OUTPUT_COLUMN_MAPPING = {
    # Basic variant info
    'chromosome': 'chromosome',
    'genomic_pos': 'genomic_pos',
    'ref': 'ref',
    'alt': 'alt',
    # Gene info
    'Symbol': 'symbol',
    'transcript': 'transcript',
    'strand': 'strand',
    # dORF overview
    'existing_dORF_count': 'existing_dORF_count',
    'affected_dORF_id': 'affected_dORF_id',
    # Classification
    'variant_type': 'variant_type',
    'dORF_type': 'dORF_type',
    # dORF coordinates
    'dORF_genomic_start': 'dORF_genomic_start',
    'dORF_genomic_end': 'dORF_genomic_end',
    'dORF_cdna_start': 'dORF_cdna_start',
    'dORF_cdna_end': 'dORF_cdna_end',
    'dORF_length': 'dORF_length',
    # Distance info
    'distance_to_CDS_stop': 'distance_to_CDS_stop',
    'overlap_CDS': 'overlap_CDS',
    # Sequence info
    'dORF_start_codon': 'dORF_start_codon',
    'dORF_stop_codon': 'dORF_stop_codon',
    'dORF_sequence': 'dORF_sequence',
    'dORF_AA': 'dORF_AA',
    'dORF_GC_percent': 'dORF_GC_percent',
    # Kozak analysis
    'kozak_sequence': 'kozak_sequence',
    'kozak_strength': 'kozak_strength',
    # Variant details
    'aa_change': 'aa_change',
    'alt_stop_exists': 'alt_stop_exists',
    'alt_stop_distance': 'alt_stop_distance',
    'indel_type': 'indel_type',
}

# Ordered column names list
OUTPUT_COLUMNS = list(OUTPUT_COLUMN_MAPPING.keys())

# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class TranscriptInfo:
    """Stores transcript structural information."""
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str = "" 
    gene_name: str = ""
    five_UTRs: List[Tuple[int, int]] = field(default_factory=list)
    CDS: List[Tuple[int, int]] = field(default_factory=list)
    three_UTRs: List[Tuple[int, int]] = field(default_factory=list)
    
    @property
    def CDS_length(self) -> int:
        return sum(end - start + 1 for start, end in self.CDS)
    
    @property
    def three_UTR_length(self) -> int:
        return sum(end - start + 1 for start, end in self.three_UTRs)
    
    @property
    def five_UTR_length(self) -> int:
        return sum(end - start + 1 for start, end in self.five_UTRs)
    
    def get_cds_stop_position(self) -> int:
        """Get the genomic position of CDS stop codon."""
        if not self.CDS:
            return 0
        if self.strand == '+':
            return max(end for _, end in self.CDS)
        else:
            return min(start for start, _ in self.CDS)


@dataclass
class dORFEntry:
    """Stores known dORF information from BED file."""
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    dorf_id: str = ""


@dataclass
class VariantEffect:
    """Stores the annotated effect of a variant on dORF."""
    # Basic info
    chromosome: str
    genomic_pos: int
    ref: str
    alt: str
    symbol: str
    transcript: str
    strand: str
    
    # dORF info
    existing_dORF_count: int = 0
    affected_dORF_id: str = ""
    variant_type: str = ""
    dORF_type: str = ""  # 'predicted' or 'existing'
    
    # Coordinates
    dORF_genomic_start: int = 0
    dORF_genomic_end: int = 0
    dORF_cdna_start: int = 0
    dORF_cdna_end: int = 0
    dORF_length: int = 0
    
    # Distance
    distance_to_CDS_stop: int = 0
    overlap_CDS: str = "No"
    
    # Sequence
    dORF_start_codon: str = ""
    dORF_stop_codon: str = ""
    dORF_sequence: str = ""
    dORF_AA: str = ""
    dORF_GC_percent: float = 0.0
    
    # Kozak
    kozak_sequence: str = ""
    kozak_strength: str = ""
    
    # Variant details
    aa_change: str = ""
    alt_stop_exists: str = "NA"
    alt_stop_distance: str = "NA"
    indel_type: str = ""
    
    def to_list(self) -> List:
        """
        Convert to list for DataFrame output.
        
        Uses OUTPUT_COLUMN_MAPPING to ensure column order matches OUTPUT_COLUMNS,
        preventing misalignment when adding/removing fields.
        """
        result = []
        for col_name in OUTPUT_COLUMNS:
            attr_name = OUTPUT_COLUMN_MAPPING.get(col_name)
            if attr_name:
                value = getattr(self, attr_name, '')
                # Special formatting for GC percent
                if attr_name == 'dORF_GC_percent' and isinstance(value, float):
                    value = f"{value:.4f}"
                result.append(value)
            else:
                result.append('')
        return result


# =============================================================================
# Configuration Class
# =============================================================================

class Config:
    """Global configuration settings."""
    def __init__(self):
        self.start_codons: Set[str] = DEFAULT_START_CODONS.copy()
        self.stop_codons: Set[str] = DEFAULT_STOP_CODONS.copy()
        self.min_orf_length: int = MIN_ORF_LENGTH
        self.max_orf_length: int = MAX_ORF_LENGTH
        self.threads: int = mp.cpu_count()
        self.chr_prefix: Optional[str] = None  # 'chr' or '' or None for auto-detect


# Global config instance
config = Config()


# =============================================================================
# Sequence Cache for Performance Optimization
# =============================================================================

class SequenceCache:
    """
    LRU cache for sequence extraction to avoid redundant FASTA reads.
    
    When processing multiple variants in the same transcript, sequences are
    often read multiple times. This cache stores recently accessed sequences
    to reduce I/O overhead.
    """
    
    def __init__(self, maxsize: int = 1000):
        """
        Initialize sequence cache.
        
        Args:
            maxsize: Maximum number of cached sequences (default 1000)
        """
        self.maxsize = maxsize
        self._cache: OrderedDict[str, str] = OrderedDict()
        self._hits = 0
        self._misses = 0
    
    def _make_key(self, transcript_id: str, seq_type: str, length: int = 0) -> str:
        """Create cache key from transcript ID and sequence type."""
        return f"{transcript_id}:{seq_type}:{length}"
    
    def get(self, transcript_id: str, seq_type: str, length: int = 0) -> Optional[str]:
        """
        Get cached sequence if available.
        
        Args:
            transcript_id: Transcript identifier
            seq_type: Type of sequence ('3utr', 'cds_tail', etc.)
            length: For variable-length sequences like cds_tail
            
        Returns:
            Cached sequence or None if not in cache
        """
        key = self._make_key(transcript_id, seq_type, length)
        if key in self._cache:
            # Move to end (most recently used)
            self._cache.move_to_end(key)
            self._hits += 1
            return self._cache[key]
        self._misses += 1
        return None
    
    def put(self, transcript_id: str, seq_type: str, sequence: str, length: int = 0) -> None:
        """
        Store sequence in cache.
        
        Args:
            transcript_id: Transcript identifier
            seq_type: Type of sequence
            sequence: The sequence to cache
            length: For variable-length sequences
        """
        key = self._make_key(transcript_id, seq_type, length)
        
        if key in self._cache:
            self._cache.move_to_end(key)
            self._cache[key] = sequence
        else:
            if len(self._cache) >= self.maxsize:
                # Remove oldest item
                self._cache.popitem(last=False)
            self._cache[key] = sequence
    
    def clear(self) -> None:
        """Clear the cache."""
        self._cache.clear()
        self._hits = 0
        self._misses = 0
    
    def stats(self) -> Dict[str, Any]:
        """Return cache statistics."""
        total = self._hits + self._misses
        hit_rate = (self._hits / total * 100) if total > 0 else 0
        return {
            'size': len(self._cache),
            'hits': self._hits,
            'misses': self._misses,
            'hit_rate': f"{hit_rate:.1f}%"
        }


# Global sequence cache instance
_sequence_cache = SequenceCache(maxsize=2000)


# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(verbose: bool = False) -> logging.Logger:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger('dORFannotator')


logger = setup_logging()


# =============================================================================
# Input File Validation
# =============================================================================

def validate_input_files(genome_file: str, vcf_file: str, gtf_file: str, 
                         dorf_file: str) -> bool:
    """Validate all input files exist and have required indices."""
    errors = []
    
    # Check genome FASTA and index
    if not os.path.isfile(genome_file):
        errors.append(f"Genome FASTA not found: {genome_file}")
    else:
        fai_file = genome_file + ".fai"
        if not os.path.isfile(fai_file):
            errors.append(f"Genome FASTA index (.fai) not found: {fai_file}")
    
    # Check VCF and index
    if not os.path.isfile(vcf_file):
        errors.append(f"VCF file not found: {vcf_file}")
    else:
        # Check for tabix index
        tbi_file = vcf_file + ".tbi"
        csi_file = vcf_file + ".csi"
        if not os.path.isfile(tbi_file) and not os.path.isfile(csi_file):
            errors.append(f"VCF index (.tbi or .csi) not found for: {vcf_file}")
    
    # Check GTF/GFF file
    if not os.path.isfile(gtf_file):
        errors.append(f"GTF/GFF file not found: {gtf_file}")
    
    # Check dORF BED file
    if not os.path.isfile(dorf_file):
        errors.append(f"dORF BED file not found: {dorf_file}")
    
    if errors:
        for err in errors:
            logger.error(err)
        return False
    
    logger.info("All input files validated successfully")
    return True


# =============================================================================
# Chromosome Naming Normalization and Region Filtering
# =============================================================================

@dataclass
class GenomicRegion:
    """Represents a genomic region for filtering."""
    chrom: str
    start: Optional[int] = None
    end: Optional[int] = None
    
    def contains(self, chrom: str, pos: int) -> bool:
        """Check if a position is within this region."""
        if not self._chrom_matches(chrom):
            return False
        if self.start is not None and pos < self.start:
            return False
        if self.end is not None and pos > self.end:
            return False
        return True
    
    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        """Check if a range overlaps with this region."""
        if not self._chrom_matches(chrom):
            return False
        if self.start is not None and end < self.start:
            return False
        if self.end is not None and start > self.end:
            return False
        return True
    
    def _chrom_matches(self, chrom: str) -> bool:
        """Check if chromosome matches (handles chr prefix differences)."""
        return chromosomes_match(self.chrom, chrom)


def parse_region_string(region_str: str) -> GenomicRegion:
    """
    Parse a region string like 'chr1:1000000-2000000' or 'chr1'.
    
    Formats:
        - 'chr1' or '1' - entire chromosome
        - 'chr1:1000000-2000000' - specific region
        - 'chr1:1000000-' - from position to end
        - 'chr1:-2000000' - from start to position
    """
    if ':' not in region_str:
        # Just chromosome
        return GenomicRegion(chrom=region_str)
    
    chrom, coords = region_str.split(':', 1)
    
    if '-' not in coords:
        # Single position? Treat as start
        pos = int(coords)
        return GenomicRegion(chrom=chrom, start=pos, end=pos)
    
    start_str, end_str = coords.split('-', 1)
    start = int(start_str) if start_str else None
    end = int(end_str) if end_str else None
    
    return GenomicRegion(chrom=chrom, start=start, end=end)


def parse_chromosome_filter(chr_arg: Optional[str]) -> Optional[Set[str]]:
    """
    Parse chromosome filter argument.
    
    Args:
        chr_arg: Comma-separated chromosome list (e.g., 'chr1,chr2' or '1,2,X')
    
    Returns:
        Set of normalized chromosome names, or None if no filter
    """
    if not chr_arg:
        return None
    
    chromosomes = set()
    for c in chr_arg.split(','):
        c = c.strip()
        if c:
            # Store both with and without chr prefix for matching
            chromosomes.add(c)
            if c.startswith('chr'):
                chromosomes.add(c[3:])
            else:
                chromosomes.add('chr' + c)
    
    return chromosomes if chromosomes else None


def normalize_chromosome(chrom: str, target_has_chr: bool) -> str:
    """Normalize chromosome naming (with or without 'chr' prefix)."""
    has_chr = chrom.startswith('chr')
    
    if target_has_chr and not has_chr:
        return 'chr' + chrom
    elif not target_has_chr and has_chr:
        return chrom[3:]
    return chrom


def get_chromosome_base(chrom: str) -> str:
    """
    Get base chromosome name without 'chr' prefix.
    
    Args:
        chrom: Chromosome name (e.g., 'chr1', '1', 'chrX')
        
    Returns:
        Base chromosome name without prefix (e.g., '1', 'X')
    """
    return chrom[3:] if chrom.startswith('chr') else chrom


def chromosomes_match(chrom1: str, chrom2: str) -> bool:
    """
    Check if two chromosome names refer to the same chromosome.
    
    Handles 'chr' prefix differences automatically.
    
    Args:
        chrom1: First chromosome name
        chrom2: Second chromosome name
        
    Returns:
        True if chromosomes match
        
    Examples:
        >>> chromosomes_match('chr1', '1')
        True
        >>> chromosomes_match('chrX', 'X')
        True
        >>> chromosomes_match('chr1', 'chr2')
        False
    """
    return get_chromosome_base(chrom1) == get_chromosome_base(chrom2)


def chromosome_in_filter(chrom: str, chr_filter: Optional[Set[str]]) -> bool:
    """
    Check if chromosome passes the filter, handling 'chr' prefix variations.
    
    Args:
        chrom: Chromosome name to check
        chr_filter: Set of allowed chromosomes (None means all pass)
        
    Returns:
        True if chromosome should be processed
        
    Examples:
        >>> chromosome_in_filter('chr1', {'1', 'chr2'})
        True
        >>> chromosome_in_filter('chrX', {'1', '2'})
        False
        >>> chromosome_in_filter('chr1', None)
        True
    """
    if not chr_filter:
        return True
    
    chrom_base = get_chromosome_base(chrom)
    return (chrom in chr_filter or 
            chrom_base in chr_filter or 
            f'chr{chrom_base}' in chr_filter)


def detect_chr_prefix(genome_file: str) -> bool:
    """Detect if genome uses 'chr' prefix."""
    try:
        fasta = pysam.FastaFile(genome_file)
        references = fasta.references
        fasta.close()
        
        # Check first few chromosomes
        for ref in references[:5]:
            if ref.startswith('chr'):
                return True
        return False
    except Exception as e:
        logger.warning(f"Could not detect chr prefix: {e}")
        return True  # Default to chr prefix


def detect_vcf_chr_format(vcf_file: str) -> bool:
    """
    Detect if VCF file uses 'chr' prefix in chromosome names.
    
    Args:
        vcf_file: Path to VCF file
        
    Returns:
        True if VCF uses 'chr' prefix, False otherwise
    """
    try:
        vcf = VCF(vcf_file)
        # Get chromosome names from VCF header
        chromosomes = vcf.seqnames
        vcf.close()
        
        if not chromosomes:
            logger.warning("Could not detect chromosome format from VCF - no chromosomes found")
            return True  # Default to chr prefix
        
        # Check first few chromosomes
        for chrom in chromosomes[:10]:
            # Skip special chromosomes
            if chrom.startswith('GL') or chrom.startswith('KI') or '_' in chrom:
                continue
            # Check if it starts with 'chr'
            if chrom.startswith('chr'):
                return True
            # If we find a numeric or X/Y chromosome without 'chr', return False
            if chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                         '21', '22', 'X', 'Y', 'MT', 'M']:
                return False
        
        # If no clear indication, default to chr prefix
        return True
    except Exception as e:
        logger.warning(f"Could not detect VCF chr format: {e}")
        return True  # Default to chr prefix


def detect_file_chr_format(file_path: str, file_type: str = 'gtf') -> bool:
    """
    Detect if GTF/GFF/BED file uses 'chr' prefix in chromosome names.
    
    Args:
        file_path: Path to file
        file_type: Type of file ('gtf', 'gff', 'bed')
        
    Returns:
        True if file uses 'chr' prefix, False otherwise
    """
    try:
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = 'rt' if file_path.endswith('.gz') else 'r'
        
        chr_count = 0
        no_chr_count = 0
        lines_checked = 0
        max_lines = 100  # Check first 100 data lines
        
        with open_func(file_path, mode) as f:
            for line in f:
                # Skip comments
                if line.startswith('#'):
                    continue
                
                # Parse chromosome (first column)
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                
                chrom = fields[0]
                
                # Skip header lines
                if chrom.lower() in ['chr', 'chrom', 'chromosome']:
                    continue
                
                # Skip special chromosomes
                if chrom.startswith('GL') or chrom.startswith('KI') or '_' in chrom:
                    continue
                
                # Count chr vs no-chr
                if chrom.startswith('chr'):
                    chr_count += 1
                elif chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                              '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                              '21', '22', 'X', 'Y', 'MT', 'M'] or chrom.startswith(tuple('123456789')):
                    no_chr_count += 1
                
                lines_checked += 1
                if lines_checked >= max_lines:
                    break
        
        # Decide based on majority
        if chr_count > no_chr_count:
            return True
        elif no_chr_count > chr_count:
            return False
        else:
            # Default to chr prefix if unclear
            return True
            
    except Exception as e:
        logger.warning(f"Could not detect chr format from {file_type} file: {e}")
        return True  # Default to chr prefix


def standardize_chromosome_name(chrom: str, target_format: bool) -> str:
    """
    Standardize chromosome name to target format.
    
    This is the main function for chromosome name conversion throughout the pipeline.
    
    Args:
        chrom: Input chromosome name
        target_format: True if target uses 'chr' prefix, False otherwise
        
    Returns:
        Standardized chromosome name
    """
    return normalize_chromosome(chrom, target_format)


def check_chromosome_format_consistency(fasta_file: str, vcf_file: str, 
                                       gtf_file: str, bed_file: str) -> Dict[str, bool]:
    """
    Check chromosome format consistency across all input files.
    
    Args:
        fasta_file: Path to FASTA file
        vcf_file: Path to VCF file
        gtf_file: Path to GTF file
        bed_file: Path to BED file
        
    Returns:
        Dictionary mapping file type to whether it uses 'chr' prefix
    """
    logger.info("Checking chromosome naming format across all input files...")
    
    formats = {}
    
    # Check each file
    formats['fasta'] = detect_chr_prefix(fasta_file)
    formats['vcf'] = detect_vcf_chr_format(vcf_file)
    formats['gtf'] = detect_file_chr_format(gtf_file, 'gtf')
    formats['bed'] = detect_file_chr_format(bed_file, 'bed')
    
    # Report formats
    logger.info("  Detected chromosome formats:")
    for file_type, has_chr in formats.items():
        prefix_str = "WITH 'chr' prefix" if has_chr else "WITHOUT 'chr' prefix"
        logger.info(f"    {file_type.upper():6s}: {prefix_str}")
    
    # Check for inconsistencies
    unique_formats = set(formats.values())
    if len(unique_formats) > 1:
        logger.warning("  ⚠ INCONSISTENT chromosome naming detected!")
        logger.warning("  The program will automatically standardize all files to match FASTA format.")
    else:
        logger.info("  ✓ All files use consistent chromosome naming format.")
    
    return formats


# =============================================================================
# GTF/GFF Parsing
# =============================================================================

def parse_gtf_attributes_fast(attr_str: str) -> Tuple[str, str]:
    """
    Fast GTF/GFF attribute parser - only extracts transcript_id and gene_name.
    
    Returns:
        Tuple of (transcript_id, gene_name)
    """
    transcript_id = ""
    gene_name = ""
    
    # Try GTF format first: key "value";
    if 'transcript_id "' in attr_str:
        # GTF format
        for part in attr_str.split(';'):
            part = part.strip()
            if part.startswith('transcript_id "'):
                transcript_id = part[15:-1].split('.')[0]  # Remove version
            elif part.startswith('gene_name "'):
                gene_name = part[11:-1]
            elif part.startswith('gene_id "') and not gene_name:
                gene_name = part[9:-1]
            
            # Early exit if we have both
            if transcript_id and gene_name:
                break
    else:
        # GFF3 format: key=value;
        for part in attr_str.split(';'):
            part = part.strip()
            if '=' not in part:
                continue
            key, _, value = part.partition('=')
            if key == 'transcript_id' or key == 'Parent':
                transcript_id = value.split('.')[0]
            elif key == 'gene_name' or key == 'Name':
                gene_name = value
            elif key == 'gene_id' and not gene_name:
                gene_name = value
            
            if transcript_id and gene_name:
                break
    
    return transcript_id, gene_name


def parse_gtf_file(gtf_file: str, 
                   chr_filter: Optional[Set[str]] = None,
                   region_filter: Optional[GenomicRegion] = None,
                   target_chr_format: Optional[bool] = None) -> Dict[str, TranscriptInfo]:
    """
    Parse GTF/GFF file to extract transcript structures.
    
    Optimized for speed:
    - Only processes relevant feature types (transcript, CDS, UTR)
    - Fast attribute parsing
    - Progress reporting for large files
    - Optional chromosome/region filtering
    - Chromosome name standardization to match FASTA format
    
    Args:
        gtf_file: Path to GTF/GFF file
        chr_filter: Set of chromosomes to include (None = all)
        region_filter: Genomic region to filter (None = all)
        target_chr_format: Target chromosome format (True=chr prefix, False=no prefix, None=auto-detect)
    """
    transcripts: Dict[str, TranscriptInfo] = {}
    gene_names: Dict[str, str] = {}
    
    # Detect GTF chromosome format if target not specified
    if target_chr_format is not None:
        gtf_has_chr = detect_file_chr_format(gtf_file, 'gtf')
        needs_conversion = (gtf_has_chr != target_chr_format)
        if needs_conversion:
            logger.info(f"  GTF uses {'chr prefix' if gtf_has_chr else 'no chr prefix'}, "
                       f"converting to {'chr prefix' if target_chr_format else 'no chr prefix'} to match FASTA")
    else:
        needs_conversion = False
    
    # Features we care about (including GFF3 variants)
    RELEVANT_FEATURES = {'transcript', 'mRNA', 'CDS', 'three_prime_UTR', 'five_prime_UTR'}
    
    filter_info = []
    if chr_filter:
        filter_info.append(f"chromosomes: {','.join(sorted(chr_filter)[:5])}...")
    if region_filter:
        region_str = f"{region_filter.chrom}"
        if region_filter.start or region_filter.end:
            region_str += f":{region_filter.start or ''}-{region_filter.end or ''}"
        filter_info.append(f"region: {region_str}")
    
    filter_msg = f" (filtering: {'; '.join(filter_info)})" if filter_info else ""
    logger.info(f"Parsing GTF/GFF file: {gtf_file}{filter_msg}")
    
    # Get file size for progress reporting
    file_size = os.path.getsize(gtf_file)
    is_gzipped = gtf_file.endswith('.gz')
    
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    bytes_read = 0
    lines_processed = 0
    last_progress = 0
    last_line_report = 0
    
    # For gzip files, estimate uncompressed size (typically 5-10x compressed size)
    # We use line-based progress reporting for gzip files
    GZIP_COMPRESSION_RATIO = 8  # Approximate ratio for text files
    estimated_size = file_size * GZIP_COMPRESSION_RATIO if is_gzipped else file_size
    
    with open_func(gtf_file, mode) as f:
        for line in f:
            if is_gzipped:
                bytes_read += len(line.encode('utf-8'))
            else:
                bytes_read += len(line)
            lines_processed += 1
            
            # Progress reporting
            if is_gzipped:
                # For gzip files, report every 500,000 lines to avoid misleading percentage
                if lines_processed >= last_line_report + 500000:
                    last_line_report = lines_processed
                    logger.info(f"  GTF parsing: {lines_processed:,} lines processed, {len(transcripts):,} transcripts found")
            else:
                # For uncompressed files, use percentage-based progress
                if file_size > 0:
                    progress = int((bytes_read / file_size) * 100)
                    if progress >= last_progress + 5:
                        last_progress = progress
                        logger.info(f"  GTF parsing progress: {progress}% ({lines_processed:,} lines, {len(transcripts):,} transcripts)")
            
            # Skip comments
            if line.startswith('#'):
                continue
            
            # Quick check for relevant features before full parsing
            # This avoids splitting lines we don't need
            has_relevant = False
            for feat in RELEVANT_FEATURES:
                if f'\t{feat}\t' in line:
                    has_relevant = True
                    break
            
            if not has_relevant:
                continue
            
            # Now parse the line
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type not in RELEVANT_FEATURES:
                continue
            
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            
            # Standardize chromosome name to target format
            if needs_conversion:
                chrom = standardize_chromosome_name(chrom, target_chr_format)
            
            # Apply chromosome filter
            if not chromosome_in_filter(chrom, chr_filter):
                continue
            
            # Apply region filter
            if region_filter and not region_filter.overlaps(chrom, start, end):
                continue
            
            # Fast attribute parsing
            transcript_id, gene_name = parse_gtf_attributes_fast(fields[8])
            
            if not transcript_id:
                continue
            
            # Store gene name
            if gene_name and transcript_id not in gene_names:
                gene_names[transcript_id] = gene_name
            
            # Initialize or update transcript
            if transcript_id not in transcripts:
                transcripts[transcript_id] = TranscriptInfo(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    transcript_id=transcript_id,
                    gene_name=gene_names.get(transcript_id, '')
                )
            
            # Add features to transcript
            tinfo = transcripts[transcript_id]
            if feature_type == 'three_prime_UTR':
                tinfo.three_UTRs.append((start, end))
            elif feature_type == 'five_prime_UTR':
                tinfo.five_UTRs.append((start, end))
            elif feature_type == 'CDS':
                tinfo.CDS.append((start, end))
    
    logger.info(f"  GTF parsing: 100% complete")
    
    # Update gene names and sort regions
    logger.info(f"  Finalizing {len(transcripts):,} transcripts...")
    
    for tid, tinfo in transcripts.items():
        if not tinfo.gene_name and tid in gene_names:
            tinfo.gene_name = gene_names[tid]
        # Sort regions in place
        tinfo.three_UTRs.sort()
        tinfo.five_UTRs.sort()
        tinfo.CDS.sort()
    
    # Filter to only keep transcripts with 3'UTR
    transcripts_with_utr = {tid: t for tid, t in transcripts.items() if t.three_UTRs}
    
    logger.info(f"Loaded {len(transcripts_with_utr):,} transcripts with 3'UTR (from {len(transcripts):,} total)")
    return transcripts_with_utr


# =============================================================================
# dORF BED File Parsing
# =============================================================================

def parse_dorf_bed(bed_file: str, 
                   chr_filter: Optional[Set[str]] = None,
                   region_filter: Optional[GenomicRegion] = None,
                   target_chr_format: Optional[bool] = None) -> Dict[str, List[dORFEntry]]:
    """Parse known dORF BED file.
    
    Format: chr start end strand transcript_id [dorf_id]
    
    Args:
        bed_file: Path to BED file
        chr_filter: Set of chromosomes to include (None = all)
        region_filter: Genomic region to filter (None = all)
        target_chr_format: Target chromosome format (True=chr prefix, False=no prefix, None=auto-detect)
    """
    dORFs: Dict[str, List[dORFEntry]] = defaultdict(list)
    
    filter_msg = " (with chromosome/region filter)" if (chr_filter or region_filter) else ""
    logger.info(f"Parsing dORF BED file: {bed_file}{filter_msg}")
    
    # Detect BED chromosome format if target specified
    if target_chr_format is not None:
        bed_has_chr = detect_file_chr_format(bed_file, 'bed')
        needs_conversion = (bed_has_chr != target_chr_format)
        if needs_conversion:
            logger.info(f"  BED uses {'chr prefix' if bed_has_chr else 'no chr prefix'}, "
                       f"converting to {'chr prefix' if target_chr_format else 'no chr prefix'} to match FASTA")
    else:
        needs_conversion = False
    
    skipped_by_filter = 0
    
    with open(bed_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('Chr'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 5:
                logger.warning(f"Skipping malformed line {line_num}: {line[:50]}...")
                continue
            
            chrom = fields[0]
            bed_start = int(fields[1])  # 0-based interbase start 
            bed_end = int(fields[2])    # 0-based interbase end
            
            # Convert to 1-based inclusive coordinates (used internally)
            # BED (start, end] -> 1-based [start+1, end]
            start_1based = bed_start + 1
            end_1based = bed_end
            
            # Standardize chromosome name to target format
            if needs_conversion:
                chrom = standardize_chromosome_name(chrom, target_chr_format)
            
            # Apply chromosome filter
            if not chromosome_in_filter(chrom, chr_filter):
                skipped_by_filter += 1
                continue
            
            # Apply region filter
            if region_filter and not region_filter.overlaps(chrom, start_1based, end_1based):
                skipped_by_filter += 1
                continue
            
            strand = fields[3]
            transcript_id = fields[4].split('.')[0]  # Remove version
            dorf_id = fields[5] if len(fields) > 5 else f"{chrom}:{start_1based}-{end_1based}"
            
            entry = dORFEntry(
                chrom=chrom,
                start=start_1based,
                end=end_1based,
                strand=strand,
                transcript_id=transcript_id,
                dorf_id=dorf_id
            )
            dORFs[transcript_id].append(entry)
    
    # Sort by position
    for entries in dORFs.values():
        entries.sort(key=lambda x: x.start)
    
    total_dorfs = sum(len(v) for v in dORFs.values())
    logger.info(f"Loaded {total_dorfs} dORFs for {len(dORFs)} transcripts")
    return dict(dORFs)


# =============================================================================
# Annotation Parsing (VEP, ANNOVAR, Custom)
# =============================================================================

def parse_vep_csq(csq_string: str, csq_format: str) -> List[Dict[str, str]]:
    """
    Parse standard VEP CSQ field.
    
    Args:
        csq_string: CSQ field value (comma-separated annotations)
        csq_format: CSQ format string from VCF header (pipe-separated field names)
    
    Returns:
        List of dicts with parsed annotation fields
    """
    results = []
    
    # Parse format string to get field names
    field_names = csq_format.split('|')
    
    # Parse each annotation (comma-separated for multi-allelic/multi-transcript)
    for annotation in csq_string.split(','):
        values = annotation.split('|')
        entry = {}
        for i, name in enumerate(field_names):
            if i < len(values):
                entry[name] = values[i]
            else:
                entry[name] = ''
        results.append(entry)
    
    return results


def parse_annovar_annotation(variant) -> Tuple[str, List[str]]:
    """
    Parse ANNOVAR annotation fields.
    
    ANNOVAR typically adds fields like:
    - Gene.refGene / Gene.ensGene
    - Func.refGene / Func.ensGene (UTR3, exonic, etc.)
    - GeneDetail.refGene (transcript details)
    - AAChange.refGene
    
    Returns:
        Tuple of (gene_symbol, list of transcript_ids)
    """
    symbol = ""
    transcript_ids = []
    
    try:
        # Try different ANNOVAR field naming conventions
        for db in ['refGene', 'ensGene', 'knownGene']:
            # Gene symbol
            gene_field = f'Gene.{db}'
            gene_value = variant.INFO.get(gene_field)
            if gene_value and not symbol:
                if isinstance(gene_value, str):
                    symbol = gene_value.split(';')[0].split('\\x3b')[0]
                else:
                    symbol = str(gene_value)
            
            # Function annotation (check if UTR3)
            func_field = f'Func.{db}'
            func_value = variant.INFO.get(func_field)
            
            # Gene detail contains transcript info
            detail_field = f'GeneDetail.{db}'
            detail_value = variant.INFO.get(detail_field)
            if detail_value:
                if isinstance(detail_value, str):
                    # Parse transcript IDs from detail (format: NM_001:c.123, NM_002:c.456)
                    for item in detail_value.replace('\\x3b', ';').split(';'):
                        if ':' in item:
                            tid = item.split(':')[0].split('.')[0]
                            if tid and tid not in transcript_ids:
                                transcript_ids.append(tid)
                        elif item.startswith(('NM_', 'NR_', 'ENST')):
                            tid = item.split('.')[0]
                            if tid not in transcript_ids:
                                transcript_ids.append(tid)
            
            # AAChange field also contains transcript info
            aachange_field = f'AAChange.{db}'
            aachange_value = variant.INFO.get(aachange_field)
            if aachange_value:
                if isinstance(aachange_value, str):
                    for item in aachange_value.replace('\\x3b', ';').split(','):
                        parts = item.split(':')
                        if len(parts) >= 2:
                            tid = parts[1].split('.')[0]
                            if tid and tid not in transcript_ids:
                                transcript_ids.append(tid)
    except KeyError as e:
        logger.debug(f"ANNOVAR field not found: {e}")
    except (ValueError, TypeError, AttributeError) as e:
        logger.debug(f"Error parsing ANNOVAR annotation: {e}")
    
    return symbol, transcript_ids


def extract_annotation_from_variant(variant, transcripts: Dict, csq_format: Optional[str] = None) -> Tuple[str, List[str]]:
    """
    Extract gene symbol and transcript IDs from variant annotation.
    
    Priority order:
    1. Custom fields (VEP_Feature, VEP_SYMBOL, VEP_Consequence)
    2. Standard VEP CSQ format
    3. ANNOVAR annotation format
    
    Args:
        variant: cyvcf2 Variant object
        transcripts: Dict of transcript_id -> TranscriptInfo
        csq_format: Optional CSQ format string from VCF header
    
    Returns:
        Tuple of (gene_symbol, list of transcript_ids)
    """
    symbol = ""
    transcript_ids = []
    
    # --- Priority 1: Custom VEP fields ---
    try:
        vep_feature = variant.INFO.get('VEP_Feature')
        vep_symbol = variant.INFO.get('VEP_SYMBOL')
        vep_consequence = variant.INFO.get('VEP_Consequence')
        
        if vep_feature:
            if isinstance(vep_feature, str):
                features = vep_feature.split('|')
            else:
                features = [str(vep_feature)]
            
            for feat in features:
                tid = feat.split('.')[0]
                if tid in transcripts and transcripts[tid].three_UTR_length > 0:
                    transcript_ids.append(tid)
            
            if vep_symbol:
                if isinstance(vep_symbol, str):
                    symbol = vep_symbol.split('|')[0]
                else:
                    symbol = str(vep_symbol)
            
            if transcript_ids:
                return symbol, transcript_ids
    except KeyError as e:
        logger.debug(f"VEP custom field not found: {e}")
    except (ValueError, TypeError) as e:
        logger.warning(f"Error parsing VEP custom fields: {e}")
    
    # --- Priority 2: Standard VEP CSQ ---
    try:
        csq_value = variant.INFO.get('CSQ')
        if csq_value and csq_format:
            if isinstance(csq_value, str):
                annotations = parse_vep_csq(csq_value, csq_format)
            else:
                annotations = parse_vep_csq(str(csq_value), csq_format)
            
            for ann in annotations:
                # Extract transcript ID (common field names: Feature, Feature_ID, Transcript)
                tid = ''
                for field in ['Feature', 'Feature_ID', 'Transcript']:
                    if field in ann and ann[field]:
                        tid = ann[field].split('.')[0]
                        break
                
                if tid and tid in transcripts and transcripts[tid].three_UTR_length > 0:
                    if tid not in transcript_ids:
                        transcript_ids.append(tid)
                
                # Extract gene symbol
                if not symbol:
                    for field in ['SYMBOL', 'Gene', 'GENE', 'Gene_Name']:
                        if field in ann and ann[field]:
                            symbol = ann[field]
                            break
            
            if transcript_ids:
                return symbol, transcript_ids
    except KeyError as e:
        logger.debug(f"CSQ field not found: {e}")
    except (ValueError, TypeError) as e:
        logger.warning(f"Error parsing VEP CSQ: {e}")
    
    # --- Priority 3: ANNOVAR format ---
    try:
        annovar_symbol, annovar_tids = parse_annovar_annotation(variant)
        if annovar_tids:
            for tid in annovar_tids:
                if tid in transcripts and transcripts[tid].three_UTR_length > 0:
                    if tid not in transcript_ids:
                        transcript_ids.append(tid)
            if not symbol and annovar_symbol:
                symbol = annovar_symbol
            
            if transcript_ids:
                return symbol, transcript_ids
    except KeyError as e:
        logger.debug(f"ANNOVAR field not found: {e}")
    except (ValueError, TypeError) as e:
        logger.warning(f"Error parsing ANNOVAR annotation: {e}")
    
    return symbol, transcript_ids


def detect_csq_format(vcf_header: str) -> Optional[str]:
    """
    Detect CSQ format string from VCF header.
    
    Looks for line like:
    ##INFO=<ID=CSQ,...,Description="...Format: Allele|Consequence|IMPACT|...">
    
    Returns:
        Format string (e.g., "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature")
        or None if not found
    """
    import re
    
    # Pattern to match CSQ format in description
    pattern = r'##INFO=<ID=CSQ,.*Description="[^"]*Format:\s*([^"]+)"'
    match = re.search(pattern, vcf_header)
    
    if match:
        return match.group(1).strip()
    
    # Alternative pattern without "Format:" prefix
    pattern2 = r'##INFO=<ID=CSQ,.*Description="([^"]+)"'
    match2 = re.search(pattern2, vcf_header)
    if match2:
        desc = match2.group(1)
        if '|' in desc:
            # Extract the pipe-separated part
            for part in desc.split(':'):
                if '|' in part:
                    return part.strip()
    
    return None


def detect_vcf_annotation_type(vcf_header: str, first_variant=None) -> str:
    """
    Detect what type of annotation the VCF file has.
    
    Checks for:
    1. Custom VEP fields (VEP_Feature, VEP_SYMBOL, VEP_Consequence)
    2. Standard VEP CSQ field
    3. ANNOVAR annotation fields
    
    Args:
        vcf_header: VCF header string
        first_variant: Optional first variant to check INFO fields
    
    Returns:
        'custom_vep', 'standard_vep', 'annovar', or 'none'
    """
    import re
    
    # Check for custom VEP fields in header
    if re.search(r'##INFO=<ID=VEP_Feature', vcf_header):
        return 'custom_vep'
    if re.search(r'##INFO=<ID=VEP_SYMBOL', vcf_header):
        return 'custom_vep'
    if re.search(r'##INFO=<ID=VEP_Consequence', vcf_header):
        return 'custom_vep'
    
    # Check for standard VEP CSQ field
    if re.search(r'##INFO=<ID=CSQ', vcf_header):
        return 'standard_vep'
    
    # Check for ANNOVAR fields
    annovar_patterns = [
        r'##INFO=<ID=Gene\.refGene',
        r'##INFO=<ID=Gene\.ensGene',
        r'##INFO=<ID=Func\.refGene',
        r'##INFO=<ID=Func\.ensGene',
        r'##INFO=<ID=AAChange\.refGene',
        r'##INFO=<ID=AAChange\.ensGene',
    ]
    for pattern in annovar_patterns:
        if re.search(pattern, vcf_header):
            return 'annovar'
    
    # Additional check: if first_variant provided, check INFO fields
    if first_variant:
        info_keys = list(first_variant.INFO.keys()) if hasattr(first_variant, 'INFO') else []
        
        # Check for custom VEP in INFO
        if any(k.startswith('VEP_') for k in info_keys):
            return 'custom_vep'
        
        # Check for CSQ in INFO
        if 'CSQ' in info_keys:
            return 'standard_vep'
        
        # Check for ANNOVAR in INFO
        if any(k.endswith('.refGene') or k.endswith('.ensGene') for k in info_keys):
            return 'annovar'
    
    return 'none'


# =============================================================================
# Sequence Utilities
# =============================================================================

def get_ordered_regions(regions: List[Tuple[int, int]], strand: str) -> List[Tuple[int, int]]:
    """Order regions in transcript 5'->3' direction."""
    if strand == '-':
        return sorted(regions, key=lambda x: x[0], reverse=True)
    return sorted(regions, key=lambda x: x[0])


def intersect_regions(region: Tuple[int, int], target_regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Calculate intersection of a region with a list of target regions.
    
    Args:
        region: (start, end) tuple, 1-based inclusive
        target_regions: List of (start, end) tuples, 1-based inclusive
    
    Returns:
        List of intersection (start, end) tuples, sorted by start position
    """
    intersections = []
    reg_start, reg_end = region
    
    for tgt_start, tgt_end in target_regions:
        # Check for overlap
        if reg_start <= tgt_end and reg_end >= tgt_start:
            # Calculate intersection
            int_start = max(reg_start, tgt_start)
            int_end = min(reg_end, tgt_end)
            if int_start <= int_end:
                intersections.append((int_start, int_end))
    
    # Sort by start position
    return sorted(intersections, key=lambda x: x[0])


def extract_sequence_from_regions(fasta: pysam.FastaFile, chrom: str, 
                                  regions: List[Tuple[int, int]], strand: str) -> str:
    """Extract and concatenate sequence from multiple regions.
    
    For negative strand: extract positive strand sequence first,
    then reverse complement the entire concatenated sequence.
    """
    if not regions:
        return ""
    
    # Always extract in genomic order first
    ordered = sorted(regions, key=lambda x: x[0])
    
    sequence = ""
    for start, end in ordered:
        try:
            # pysam uses 0-based, half-open coordinates
            seq = fasta.fetch(chrom, start - 1, end)
            sequence += seq.upper()
        except Exception as e:
            logger.warning(f"Could not fetch sequence {chrom}:{start}-{end}: {e}")
            return ""
    
    # Reverse complement for negative strand
    if strand == '-':
        sequence = str(reverse_complement(sequence))
    
    return sequence


def get_3utr_sequence(transcript: TranscriptInfo, fasta: pysam.FastaFile, 
                      use_cache: bool = True) -> str:
    """
    Get the 3'UTR sequence in transcript orientation.
    
    Args:
        transcript: Transcript information
        fasta: Open pysam FastaFile
        use_cache: Whether to use sequence cache (default True)
        
    Returns:
        3'UTR sequence in transcript orientation
    """
    # Check cache first
    if use_cache:
        cached = _sequence_cache.get(transcript.transcript_id, '3utr')
        if cached is not None:
            return cached
    
    # Extract sequence
    sequence = extract_sequence_from_regions(
        fasta, transcript.chrom, transcript.three_UTRs, transcript.strand
    )
    
    # Store in cache
    if use_cache and sequence:
        _sequence_cache.put(transcript.transcript_id, '3utr', sequence)
    
    return sequence


# =============================================================================
# Kozak Sequence Analysis
# =============================================================================

def classify_kozak_strength(kozak_seq: str, start_codon: str) -> str:
    """
    Classify Kozak strength based on Marilyn Kozak rules.
    
    - Strong: -3 is A/G AND +4 is G
    - Moderate: -3 is A/G OR +4 is G
    - Weak: neither condition met
    - '': not a valid start codon
    """
    if len(kozak_seq) < 7:
        return ''
    
    # Only analyze for ATG (canonical start)
    if start_codon != 'ATG':
        return ''
    
    minus_3 = kozak_seq[0].upper()  # Position -3
    plus_4 = kozak_seq[6].upper()   # Position +4 (A of ATG is +1)
    
    if minus_3 == 'N' or plus_4 == 'N':
        return ''
    
    is_minus3_ag = minus_3 in ('A', 'G')
    is_plus4_g = plus_4 == 'G'
    
    if is_minus3_ag and is_plus4_g:
        return 'Strong'
    elif is_minus3_ag or is_plus4_g:
        return 'Moderate'
    else:
        return 'Weak'


def get_cds_tail_sequence(transcript: TranscriptInfo, fasta: pysam.FastaFile, 
                          length: int = 10, use_cache: bool = True) -> str:
    """
    Get the last `length` bases of CDS (immediately before 3'UTR).
    
    For Kozak sequence analysis, we need bases upstream of 3'UTR when 
    dORF start codon is near the beginning of 3'UTR.
    
    Args:
        transcript: Transcript information with CDS regions
        fasta: Open pysam FastaFile
        length: Number of bases to fetch from CDS end
        use_cache: Whether to use sequence cache (default True)
    
    Returns:
        CDS tail sequence in transcript orientation (5'->3')
    """
    if not transcript.CDS or length <= 0:
        return ""
    
    # Check cache first
    if use_cache:
        cached = _sequence_cache.get(transcript.transcript_id, 'cds_tail', length)
        if cached is not None:
            return cached
    
    # Get ordered CDS regions (5'->3' in transcript)
    ordered_cds = get_ordered_regions(transcript.CDS, transcript.strand)
    
    # Collect sequence from the 3' end of CDS
    tail_seq = ""
    remaining = length
    
    # Iterate from 3' end of CDS (last exon first in transcript order)
    for start, end in reversed(ordered_cds):
        try:
            # Get sequence from genome (always positive strand first)
            seq = fasta.fetch(transcript.chrom, start - 1, end).upper()
            
            # Reverse complement for negative strand
            if transcript.strand == '-':
                seq = str(reverse_complement(seq))
            
            if len(seq) >= remaining:
                tail_seq = seq[-remaining:] + tail_seq
                break
            else:
                tail_seq = seq + tail_seq
                remaining -= len(seq)
        except KeyError:
            logger.debug(f"Chromosome {transcript.chrom} not found in FASTA for CDS tail")
            break
        except ValueError as e:
            logger.debug(f"Invalid CDS coordinates {start}-{end}: {e}")
            continue
        except Exception as e:
            logger.warning(f"Error fetching CDS sequence {transcript.chrom}:{start}-{end}: {e}")
            continue
    
    # Get the final sequence
    result = tail_seq[-length:] if len(tail_seq) >= length else tail_seq
    
    # Store in cache
    if use_cache and result:
        _sequence_cache.put(transcript.transcript_id, 'cds_tail', result, length)
    
    return result


def get_kozak_context(utr_seq: str, start_idx: int, cds_tail: str = "") -> Tuple[str, str]:
    """
    Get 7bp Kozak context (-3 to +4, where A of ATG is +1).
    
    Args:
        utr_seq: 3'UTR sequence
        start_idx: 0-based index of start codon in UTR
        cds_tail: Last few bases of CDS (for positions before UTR start)
    
    Returns:
        Tuple of (kozak_sequence, kozak_strength)
    
    Position layout:
        -3  -2  -1  +1  +2  +3  +4
        [CDS tail ] [A   T   G] [next base]
                    [start codon]
    
    If start_idx < 3, we need bases from CDS tail for positions -3, -2, -1
    """
    kozak_bases = []
    
    for offset in range(-3, 4):
        pos = start_idx + offset
        if 0 <= pos < len(utr_seq):
            # Position is within UTR
            kozak_bases.append(utr_seq[pos])
        elif pos < 0:
            # Position is before UTR start, need CDS tail
            if cds_tail:
                # cds_tail is the last N bases of CDS
                # pos=-1 means we need the last base of CDS (cds_tail[-1])
                # pos=-2 means we need the second to last (cds_tail[-2])
                # pos=-3 means we need the third to last (cds_tail[-3])
                cds_idx = len(cds_tail) + pos  # e.g., if len=10 and pos=-1, idx=9
                if 0 <= cds_idx < len(cds_tail):
                    kozak_bases.append(cds_tail[cds_idx])
                else:
                    kozak_bases.append('N')
            else:
                kozak_bases.append('N')
        else:
            # Position beyond UTR end
            kozak_bases.append('N')
    
    kozak_seq = ''.join(kozak_bases).upper()
    start_codon = utr_seq[start_idx:start_idx+3].upper() if start_idx + 3 <= len(utr_seq) else ''
    
    return kozak_seq, classify_kozak_strength(kozak_seq, start_codon)


# =============================================================================
# ORF Finding and Validation
# =============================================================================

def find_orf_from_start(sequence: str, start_idx: int) -> str:
    """Find complete ORF from a start codon position."""
    if start_idx < 0 or start_idx + 3 > len(sequence):
        return ""
    
    if sequence[start_idx:start_idx+3] not in config.start_codons:
        return ""
    
    orf = sequence[start_idx:start_idx+3]
    pos = start_idx + 3
    
    while pos + 3 <= len(sequence):
        codon = sequence[pos:pos+3]
        orf += codon
        if codon in config.stop_codons:
            return orf
        pos += 3
    
    return ""  # No stop codon found


def find_upstream_start(sequence: str, stop_idx: int) -> Tuple[str, int]:
    """
    Find the most upstream (5') start codon in the same reading frame.
    
    Returns:
        Tuple of (orf_sequence, start_index) or ("", -1) if not found
    """
    if stop_idx < 0 or stop_idx + 3 > len(sequence):
        return "", -1
    
    # Find most upstream start in same frame
    best_start = -1
    frame = stop_idx % 3
    
    for pos in range(frame, stop_idx, 3):
        if pos + 3 <= len(sequence):
            codon = sequence[pos:pos+3]
            if codon in config.start_codons:
                # Check no premature stop between this start and the final stop
                has_premature_stop = False
                for check_pos in range(pos + 3, stop_idx, 3):
                    if sequence[check_pos:check_pos+3] in config.stop_codons:
                        has_premature_stop = True
                        break
                
                if not has_premature_stop:
                    best_start = pos
                    break  # Found most upstream valid start
    
    if best_start >= 0:
        orf_seq = sequence[best_start:stop_idx+3]
        return orf_seq, best_start
    
    return "", -1


def validate_orf(orf_seq: str) -> bool:
    """Validate ORF meets length and structure criteria."""
    if not orf_seq:
        return False
    
    length = len(orf_seq)
    
    return (
        length >= config.min_orf_length and
        length <= config.max_orf_length and
        length % 3 == 0 and
        orf_seq[:3] in config.start_codons and
        orf_seq[-3:] in config.stop_codons
    )


def find_alternative_stop(sequence: str, start_after: int) -> Tuple[int, int]:
    """
    Find alternative stop codon downstream.
    
    Returns:
        Tuple of (stop_index, distance) or (-1, -1) if not found
    """
    frame = start_after % 3
    
    for pos in range(start_after + 3, len(sequence) - 2, 3):
        if pos % 3 == frame:
            if sequence[pos:pos+3] in config.stop_codons:
                return pos, pos - start_after
    
    return -1, -1


# =============================================================================
# Variant Application
# =============================================================================

def apply_variant_to_sequence(sequence: str, pos_0based: int, ref: str, alt: str) -> Optional[str]:
    """
    Apply a variant to a sequence.
    
    Args:
        sequence: Original sequence
        pos_0based: 0-based position in sequence
        ref: Reference allele
        alt: Alternate allele
    
    Returns:
        Modified sequence or None if invalid
    """
    if pos_0based < 0 or pos_0based + len(ref) > len(sequence):
        return None
    
    # Verify reference matches
    actual_ref = sequence[pos_0based:pos_0based + len(ref)]
    if actual_ref.upper() != ref.upper():
        logger.debug(f"Reference mismatch: expected {ref}, got {actual_ref}")
        # Continue anyway, as this might be due to strand issues
    
    # Apply the variant
    return sequence[:pos_0based] + alt + sequence[pos_0based + len(ref):]


# =============================================================================
# Coordinate Conversion
# =============================================================================

def genomic_to_cdna_position(transcript: TranscriptInfo, genomic_pos: int) -> int:
    """
    Convert genomic position to cDNA position within 3'UTR.
    
    Returns:
        1-based cDNA position, or 0 if not in 3'UTR
    """
    cdna_pos = 0
    
    ordered_utrs = get_ordered_regions(transcript.three_UTRs, transcript.strand)
    
    for start, end in ordered_utrs:
        if start <= genomic_pos <= end:
            if transcript.strand == '+':
                cdna_pos += (genomic_pos - start + 1)
            else:
                cdna_pos += (end - genomic_pos + 1)
            return cdna_pos
        else:
            cdna_pos += (end - start + 1)
    
    return 0  # Position not in 3'UTR


def cdna_to_genomic_position(transcript: TranscriptInfo, cdna_pos: int) -> int:
    """
    Convert cDNA position to genomic position.
    
    Args:
        cdna_pos: 1-based cDNA position
    
    Returns:
        Genomic position (1-based), or 0 if invalid
    """
    remaining = cdna_pos
    ordered_utrs = get_ordered_regions(transcript.three_UTRs, transcript.strand)
    
    for start, end in ordered_utrs:
        region_len = end - start + 1
        if remaining <= region_len:
            if transcript.strand == '+':
                return start + remaining - 1
            else:
                return end - remaining + 1
        remaining -= region_len
    
    return 0


def cdna_to_genomic_position_in_cds(transcript: TranscriptInfo, cds_offset: int) -> int:
    """
    Convert a position relative to CDS end (3' end of CDS) to genomic position.
    
    Args:
        transcript: TranscriptInfo object
        cds_offset: Number of bp upstream from the 3' end of CDS (1-based)
                   e.g., cds_offset=1 means the last bp of CDS
    
    Returns:
        Genomic position (1-based), or 0 if invalid
    """
    if not transcript.CDS or cds_offset <= 0:
        return 0
    
    # Get CDS regions ordered in 5'->3' direction (transcript orientation)
    ordered_cds = get_ordered_regions(transcript.CDS, transcript.strand)
    
    # We need to count from 3' end of CDS, so reverse the order
    remaining = cds_offset
    
    for start, end in reversed(ordered_cds):
        region_len = end - start + 1
        if remaining <= region_len:
            if transcript.strand == '+':
                # For + strand, 3' end is at high coordinate
                return end - remaining + 1
            else:
                # For - strand, 3' end is at low coordinate
                return start + remaining - 1
        remaining -= region_len
    
    return 0


# =============================================================================
# Variant Classification Functions
# =============================================================================

def analyze_dStart_gained(ref_utr: str, alt_utr: str, variant_cdna_pos: int,
                         transcript: TranscriptInfo, ref: str, alt: str) -> List[Dict]:
    """
    Check if variant creates a new start codon.
    
    Look for new start codons around the variant position, then find
    the first in-frame stop codon downstream.
    """
    results = []
    
    if not alt_utr:
        return results
    
    # Check positions around the variant for new start codons
    for offset in range(-2, len(alt) + 1):
        start_idx = variant_cdna_pos - 1 + offset
        
        if start_idx < 0 or start_idx + 3 > len(alt_utr):
            continue
        
        new_codon = alt_utr[start_idx:start_idx+3]
        
        if new_codon in config.start_codons:
            # Check if this was NOT a start codon in reference
            ref_idx = start_idx
            if ref_idx >= 0 and ref_idx + 3 <= len(ref_utr):
                old_codon = ref_utr[ref_idx:ref_idx+3]
                if old_codon in config.start_codons:
                    continue  # Not a new start
            
            # Find complete ORF
            orf_seq = find_orf_from_start(alt_utr, start_idx)
            
            if validate_orf(orf_seq):
                results.append({
                    'type': 'dStart_gained',
                    'orf_seq': orf_seq,
                    'cdna_start': start_idx + 1,
                    'cdna_end': start_idx + len(orf_seq),
                    'is_new': True
                })
    
    return results


def analyze_dStop_gained(ref_utr: str, alt_utr: str, variant_cdna_pos: int,
                        transcript: TranscriptInfo, ref: str, alt: str,
                        cds_tail_seq: str = "", cds_tail_length: int = 0) -> List[Dict]:
    """
    Check if variant creates a new stop codon.
    
    Find the most upstream start codon in the same reading frame.
    Search extends into CDS region (up to 300bp) to find longest possible dORF.
    
    Args:
        ref_utr: Reference 3'UTR sequence
        alt_utr: Alternate 3'UTR sequence with variant applied
        variant_cdna_pos: 1-based position of variant in 3'UTR cDNA coordinates
        transcript: TranscriptInfo object
        ref: Reference allele
        alt: Alternate allele
        cds_tail_seq: Last portion of CDS sequence (up to 300bp) for upstream search
        cds_tail_length: Length of CDS tail sequence provided
    """
    results = []
    
    if not alt_utr:
        return results
    
    # Combine CDS tail and UTR for extended search
    # CDS tail is in transcript orientation (5'->3'), so prepend to UTR
    extended_alt_seq = cds_tail_seq + alt_utr if cds_tail_seq else alt_utr
    extended_ref_seq = cds_tail_seq + ref_utr if cds_tail_seq else ref_utr
    
    # Offset for coordinates: positions in extended_seq need this offset subtracted
    # to get positions in alt_utr
    cds_offset = len(cds_tail_seq)
    
    # Adjust variant position for extended sequence
    extended_variant_pos = variant_cdna_pos + cds_offset
    
    # Check positions around the variant for new stop codons
    for offset in range(-2, len(alt) + 1):
        stop_idx_extended = extended_variant_pos - 1 + offset
        
        if stop_idx_extended < 0 or stop_idx_extended + 3 > len(extended_alt_seq):
            continue
        
        new_codon = extended_alt_seq[stop_idx_extended:stop_idx_extended+3]
        
        if new_codon in config.stop_codons:
            # Check if this was NOT a stop codon in reference
            ref_idx = stop_idx_extended
            if ref_idx >= 0 and ref_idx + 3 <= len(extended_ref_seq):
                old_codon = extended_ref_seq[ref_idx:ref_idx+3]
                if old_codon in config.stop_codons:
                    continue  # Not a new stop
            
            # Find upstream start in extended sequence (includes CDS)
            orf_seq, start_idx_extended = find_upstream_start_extended(
                extended_alt_seq, stop_idx_extended, cds_offset
            )
            
            if validate_orf(orf_seq):
                # Calculate coordinates in UTR space
                # start_idx_extended is 0-based position in extended sequence
                # cdna_start should be 1-based position in UTR
                # If start is in CDS (start_idx_extended < cds_offset), cdna_start will be negative
                cdna_start = start_idx_extended - cds_offset + 1
                cdna_end = stop_idx_extended - cds_offset + 3
                
                # Check if dORF overlaps with CDS
                overlap_cds = start_idx_extended < cds_offset
                
                results.append({
                    'type': 'dStop_gained',
                    'orf_seq': orf_seq,
                    'cdna_start': cdna_start,
                    'cdna_end': cdna_end,
                    'is_new': True,
                    'overlap_cds': overlap_cds,
                    'cds_overlap_length': max(0, cds_offset - start_idx_extended) if overlap_cds else 0
                })
    
    return results


def find_upstream_start_extended(sequence: str, stop_idx: int, cds_offset: int) -> Tuple[str, int]:
    """
    Find the most upstream (5') start codon in the same reading frame.
    Search is limited to produce ORF within max_orf_length (default 303bp).
    Prioritizes finding the longest valid ORF.
    
    Args:
        sequence: Extended sequence (CDS tail + 3'UTR)
        stop_idx: 0-based index of stop codon in extended sequence
        cds_offset: Length of CDS portion in extended sequence
    
    Returns:
        Tuple of (orf_sequence, start_index) or ("", -1) if not found
    """
    if stop_idx < 0 or stop_idx + 3 > len(sequence):
        return "", -1
    
    # Find most upstream start in same frame that produces valid ORF
    best_start = -1
    frame = stop_idx % 3
    
    # Calculate search limit based on max ORF length
    # Stop codon at stop_idx, so start must be at least 3 codons before (min ORF length)
    min_start = max(0, stop_idx - config.max_orf_length + 3)
    
    for pos in range(min_start + (frame - min_start % 3) % 3, stop_idx, 3):
        if pos + 3 <= len(sequence):
            codon = sequence[pos:pos+3]
            if codon in config.start_codons:
                # Check no premature stop between this start and the final stop
                has_premature_stop = False
                for check_pos in range(pos + 3, stop_idx, 3):
                    if sequence[check_pos:check_pos+3] in config.stop_codons:
                        has_premature_stop = True
                        break
                
                if not has_premature_stop:
                    orf_len = stop_idx + 3 - pos
                    if config.min_orf_length <= orf_len <= config.max_orf_length:
                        best_start = pos
                        break  # Found most upstream valid start
    
    if best_start >= 0:
        orf_seq = sequence[best_start:stop_idx+3]
        return orf_seq, best_start
    
    return "", -1


def analyze_existing_dorf_effects(ref_utr: str, alt_utr: str, variant_cdna_pos: int,
                                   dorf_entry: dORFEntry, orf_seq: str, orf_start_cdna: int,
                                   ref: str, alt: str) -> List[Dict]:
    """
    Analyze effects on an existing dORF.
    
    Returns list of effect dictionaries.
    """
    results = []
    
    if not alt_utr or not orf_seq:
        return results
    
    orf_end_cdna = orf_start_cdna + len(orf_seq) - 1
    
    # Check if variant affects this dORF
    variant_end_cdna = variant_cdna_pos + len(ref) - 1
    if variant_end_cdna < orf_start_cdna or variant_cdna_pos > orf_end_cdna:
        return results  # Variant outside this dORF
    
    # Calculate position within ORF (0-based)
    pos_in_orf = variant_cdna_pos - orf_start_cdna
    
    # Apply variant to get modified ORF
    # First, extract the ORF portion from alt_utr
    orf_start_0based = orf_start_cdna - 1
    
    # Handle length changes from indels
    len_diff = len(alt) - len(ref)
    new_orf_len = len(orf_seq) + len_diff
    
    if orf_start_0based < 0 or orf_start_0based >= len(alt_utr):
        return results
    
    # Get the modified ORF region from alt_utr
    modified_orf = alt_utr[orf_start_0based:orf_start_0based + new_orf_len]
    
    if len(modified_orf) < 3:
        return results
    
    # Check for start codon effects
    original_start = orf_seq[:3]
    modified_start = modified_orf[:3] if len(modified_orf) >= 3 else ""
    
    if pos_in_orf >= 0 and pos_in_orf < 3:  # Variant affects start codon
        if original_start in config.start_codons:
            if modified_start not in config.start_codons:
                # Start lost
                results.append({
                    'type': 'dStart_lost',
                    'orf_seq': orf_seq,
                    'cdna_start': orf_start_cdna,
                    'cdna_end': orf_end_cdna,
                    'dorf_id': dorf_entry.dorf_id
                })
            elif modified_start != original_start:
                # Start changed
                results.append({
                    'type': 'dStart_change',
                    'orf_seq': orf_seq,
                    'cdna_start': orf_start_cdna,
                    'cdna_end': orf_end_cdna,
                    'dorf_id': dorf_entry.dorf_id,
                    'old_start': original_start,
                    'new_start': modified_start
                })
        return results  # Start codon effects take priority
    
    # Check for stop codon effects
    stop_pos_in_orf = len(orf_seq) - 3
    original_stop = orf_seq[-3:]
    
    if pos_in_orf >= stop_pos_in_orf and pos_in_orf < len(orf_seq):
        modified_stop = modified_orf[-3:] if len(modified_orf) >= 3 else ""
        
        if original_stop in config.stop_codons:
            if modified_stop not in config.stop_codons:
                # Stop lost - find alternative stop
                alt_stop_idx, alt_stop_dist = find_alternative_stop(
                    alt_utr, orf_start_0based + len(orf_seq) - 3
                )
                results.append({
                    'type': 'dStop_lost',
                    'orf_seq': orf_seq,
                    'cdna_start': orf_start_cdna,
                    'cdna_end': orf_end_cdna,
                    'dorf_id': dorf_entry.dorf_id,
                    'alt_stop_exists': alt_stop_idx >= 0,
                    'alt_stop_distance': alt_stop_dist if alt_stop_idx >= 0 else None
                })
            elif modified_stop != original_stop:
                # Stop changed
                results.append({
                    'type': 'dStop_change',
                    'orf_seq': orf_seq,
                    'cdna_start': orf_start_cdna,
                    'cdna_end': orf_end_cdna,
                    'dorf_id': dorf_entry.dorf_id,
                    'old_stop': original_stop,
                    'new_stop': modified_stop
                })
        return results  # Stop codon effects take priority
    
    # Check for indel effects (frameshift vs inframe)
    if len(ref) != len(alt):
        indel_len = abs(len(alt) - len(ref))
        if indel_len % 3 == 0:
            # Inframe indel
            results.append({
                'type': 'dIndel_inframe',
                'orf_seq': orf_seq,
                'cdna_start': orf_start_cdna,
                'cdna_end': orf_end_cdna,
                'dorf_id': dorf_entry.dorf_id,
                'indel_length': indel_len,
                'is_insertion': len(alt) > len(ref)
            })
        else:
            # Frameshift
            results.append({
                'type': 'dIndel_frameshift',
                'orf_seq': orf_seq,
                'cdna_start': orf_start_cdna,
                'cdna_end': orf_end_cdna,
                'dorf_id': dorf_entry.dorf_id,
                'indel_length': indel_len,
                'is_insertion': len(alt) > len(ref)
            })
        return results
    
    # Check for missense (SNV causing amino acid change)
    if len(ref) == 1 and len(alt) == 1:
        # Translate both (ensure sequence length is multiple of 3 to avoid Biopython warning)
        try:
            ref_protein = str(translate(orf_seq)) if len(orf_seq) % 3 == 0 else ""
            alt_protein = str(translate(modified_orf)) if len(modified_orf) % 3 == 0 else ""
            
            if ref_protein and alt_protein and ref_protein != alt_protein:
                # Find changed position
                for i, (r, a) in enumerate(zip(ref_protein, alt_protein)):
                    if r != a:
                        results.append({
                            'type': 'dMissense',
                            'orf_seq': orf_seq,
                            'cdna_start': orf_start_cdna,
                            'cdna_end': orf_end_cdna,
                            'dorf_id': dorf_entry.dorf_id,
                            'aa_change': f"p.{r}{i+1}{a}"
                        })
                        break
        except Exception:
            pass
    
    return results


# =============================================================================
# Main Variant Processing
# =============================================================================

def process_variant(chrom: str, pos: int, ref: str, alt: str,
                   transcript: TranscriptInfo, fasta: pysam.FastaFile,
                   existing_dorfs: List[dORFEntry], symbol: str) -> List[VariantEffect]:
    """
    Process a single variant for a transcript.
    
    Handles negative strand properly by:
    1. Getting reference 3'UTR sequence
    2. Applying variant in genomic coordinates
    3. Getting modified 3'UTR sequence
    4. Analyzing effects
    """
    results = []
    
    # Get reference UTR sequence
    ref_utr = get_3utr_sequence(transcript, fasta)
    if not ref_utr:
        return results
    
    # Get CDS tail for Kozak sequence analysis (need at least 3 bases for -3 position)
    cds_tail_kozak = get_cds_tail_sequence(transcript, fasta, length=10)
    
    # Get extended CDS tail for dStop_gained analysis (up to 300bp for searching upstream starts)
    cds_tail_extended = get_cds_tail_sequence(transcript, fasta, length=300)
    
    # Convert genomic position to cDNA position
    cdna_pos = genomic_to_cdna_position(transcript, pos)
    if cdna_pos == 0:
        return results  # Variant not in 3'UTR
    
    # For negative strand, we need to reverse complement the alleles
    if transcript.strand == '-':
        ref_rc = str(reverse_complement(ref))
        alt_rc = str(reverse_complement(alt))
    else:
        ref_rc = ref
        alt_rc = alt
    if transcript.strand == '-' and len(ref) > 1:
        # For negative strand multi-base variants:
        # cdna_pos points to the 3' end of variant in transcript coordinates
        # Adjust to the 5' end (which is len(ref)-1 positions earlier in cDNA)
        apply_start_0based = cdna_pos - len(ref)
    else:
        # For positive strand, or single-base variants on negative strand
        apply_start_0based = cdna_pos - 1
    
    # Validate the apply position is within bounds
    if apply_start_0based < 0:
        logger.debug(f"Variant application position out of bounds: {apply_start_0based}")
        return results
    
    # Apply variant to get modified UTR
    alt_utr = apply_variant_to_sequence(ref_utr, apply_start_0based, ref_rc, alt_rc)
    
    if not alt_utr:
        return results
    
    # Calculate the 1-based cDNA position of variant's 5' end in transcript
    # This is used for the analysis functions to correctly locate codon positions
    if transcript.strand == '-' and len(ref) > 1:
        # For negative strand multi-base: 5' end is at lower cDNA position
        variant_cdna_5prime = cdna_pos - len(ref) + 1
    else:
        variant_cdna_5prime = cdna_pos
    
    # Count existing dORFs
    dorf_count = len(existing_dorfs)
    
    # Base effect template
    base_effect = VariantEffect(
        chromosome=chrom,
        genomic_pos=pos,
        ref=ref,
        alt=alt,
        symbol=symbol,
        transcript=transcript.gene_name if transcript.gene_name else "",
        strand=transcript.strand,
        existing_dORF_count=dorf_count
    )
    
    effects_found = []
    
    # 1. Check for gained effects (new dORFs)
    # Use variant_cdna_5prime for correct codon position searching
    gained_starts = analyze_dStart_gained(ref_utr, alt_utr, variant_cdna_5prime, transcript, ref_rc, alt_rc)
    gained_stops = analyze_dStop_gained(ref_utr, alt_utr, variant_cdna_5prime, transcript, ref_rc, alt_rc,
                                        cds_tail_extended, len(cds_tail_extended))
    
    for effect_data in gained_starts + gained_stops:
        effect = VariantEffect(**{k: v for k, v in base_effect.__dict__.items()})
        effect.variant_type = effect_data['type']
        effect.dORF_type = 'predicted'
        effect.dORF_cdna_start = effect_data['cdna_start']
        effect.dORF_cdna_end = effect_data['cdna_end']
        effect.dORF_length = len(effect_data['orf_seq'])
        effect.dORF_sequence = effect_data['orf_seq']
        effect.dORF_start_codon = effect_data['orf_seq'][:3]
        effect.dORF_stop_codon = effect_data['orf_seq'][-3:]
        
        # Set overlap_CDS based on analysis results
        if effect_data.get('overlap_cds', False):
            effect.overlap_CDS = "TRUE"
        else:
            effect.overlap_CDS = "FALSE"
        
        try:
            orf_seq = effect_data['orf_seq']
            # Ensure sequence length is multiple of 3 before translation
            if len(orf_seq) >= 3:
                if len(orf_seq) % 3 == 0:
                    effect.dORF_AA = str(translate(orf_seq))
                else:
                    # Trim to multiple of 3
                    trimmed_len = (len(orf_seq) // 3) * 3
                    effect.dORF_AA = str(translate(orf_seq[:trimmed_len]))
            effect.dORF_GC_percent = gc_fraction(orf_seq)
        except ValueError as e:
            logger.debug(f"Translation failed for ORF sequence: {e}")
        except Exception as e:
            logger.warning(f"Error calculating ORF properties: {e}")
        
        # Get genomic coordinates
        # For dStop_gained with CDS overlap, cdna_start can be negative (in CDS)
        if effect.dORF_cdna_start > 0:
            effect.dORF_genomic_start = cdna_to_genomic_position(transcript, effect.dORF_cdna_start)
        else:
            # Start is in CDS region
            cds_offset = 1 - effect.dORF_cdna_start  # How many bp into CDS
            effect.dORF_genomic_start = cdna_to_genomic_position_in_cds(transcript, cds_offset)
        
        effect.dORF_genomic_end = cdna_to_genomic_position(transcript, effect.dORF_cdna_end)
        
        # Calculate distance to CDS
        cds_stop = transcript.get_cds_stop_position()
        if cds_stop and effect.dORF_genomic_start:
            effect.distance_to_CDS_stop = abs(effect.dORF_genomic_start - cds_stop)
        
        # Kozak analysis (pass CDS tail for positions near UTR start)
        # For positions in CDS (negative cdna_start), use the extended CDS tail
        if effect.dORF_cdna_start <= 0:
            # Start is in CDS, need to calculate Kozak from CDS + UTR
            start_in_extended = len(cds_tail_extended) + effect.dORF_cdna_start - 1
            extended_seq = cds_tail_extended + alt_utr
            kozak_seq, kozak_strength = get_kozak_context(extended_seq, start_in_extended, "")
        else:
            kozak_seq, kozak_strength = get_kozak_context(alt_utr, effect.dORF_cdna_start - 1, cds_tail_kozak)
        effect.kozak_sequence = kozak_seq
        effect.kozak_strength = kozak_strength
        
        effects_found.append(effect)
    
    # 2. Analyze effects on existing dORFs
    for dorf_entry in existing_dorfs:
        # Get dORF sequence from reference, intersecting with 3'UTR regions
        # This ensures we only extract sequence from 3'UTR exons, not introns or other regions
        dorf_region = (dorf_entry.start, dorf_entry.end)
        dorf_utr_regions = intersect_regions(dorf_region, transcript.three_UTRs)
        
        if not dorf_utr_regions:
            # dORF doesn't overlap with any 3'UTR region, skip
            logger.debug(f"dORF {dorf_entry.dorf_id} does not overlap with 3'UTR of transcript, skipping")
            continue
        
        # Extract sequence only from 3'UTR intersections
        dorf_seq = extract_sequence_from_regions(fasta, dorf_entry.chrom, dorf_utr_regions, dorf_entry.strand)
        
        if not dorf_seq:
            continue
        
        # Update dORF entry with actual 3'UTR-intersected boundaries for coordinate calculation
        if transcript.strand == '+':
            actual_dorf_start = min(start for start, _ in dorf_utr_regions)
            actual_dorf_end = max(end for _, end in dorf_utr_regions)
        else:
            actual_dorf_start = min(start for start, _ in dorf_utr_regions)
            actual_dorf_end = max(end for _, end in dorf_utr_regions)
        
        # Calculate dORF cDNA position using the actual 3'UTR-intersected boundaries
        if transcript.strand == '+':
            dorf_cdna_start = genomic_to_cdna_position(transcript, actual_dorf_start)
        else:
            dorf_cdna_start = genomic_to_cdna_position(transcript, actual_dorf_end)
        
        existing_effects = analyze_existing_dorf_effects(
            ref_utr, alt_utr, variant_cdna_5prime, dorf_entry, dorf_seq, dorf_cdna_start, ref_rc, alt_rc
        )
        
        for effect_data in existing_effects:
            effect = VariantEffect(**{k: v for k, v in base_effect.__dict__.items()})
            effect.variant_type = effect_data['type']
            effect.dORF_type = 'existing'
            effect.affected_dORF_id = effect_data.get('dorf_id', '')
            effect.dORF_cdna_start = effect_data['cdna_start']
            effect.dORF_cdna_end = effect_data['cdna_end']
            effect.dORF_length = len(effect_data['orf_seq'])
            effect.dORF_sequence = effect_data['orf_seq']
            effect.dORF_start_codon = effect_data['orf_seq'][:3]
            effect.dORF_stop_codon = effect_data['orf_seq'][-3:]
            
            try:
                orf_seq = effect_data['orf_seq']
                # Ensure sequence length is multiple of 3 to avoid Biopython warning
                if len(orf_seq) >= 3:
                    if len(orf_seq) % 3 == 0:
                        effect.dORF_AA = str(translate(orf_seq))
                    else:
                        trimmed_len = (len(orf_seq) // 3) * 3
                        effect.dORF_AA = str(translate(orf_seq[:trimmed_len])) if trimmed_len > 0 else ""
                effect.dORF_GC_percent = gc_fraction(orf_seq)
            except ValueError as e:
                logger.debug(f"Translation failed for existing dORF: {e}")
            except Exception as e:
                logger.warning(f"Error calculating existing dORF properties: {e}")
            
            # Use the actual 3'UTR-intersected boundaries, not original BED coordinates
            effect.dORF_genomic_start = actual_dorf_start
            effect.dORF_genomic_end = actual_dorf_end
            
            # Distance to CDS
            cds_stop = transcript.get_cds_stop_position()
            if cds_stop:
                effect.distance_to_CDS_stop = abs(actual_dorf_start - cds_stop)
            
            # Kozak (pass CDS tail for positions near UTR start)
            kozak_seq, kozak_strength = get_kozak_context(ref_utr, effect.dORF_cdna_start - 1, cds_tail_kozak)
            effect.kozak_sequence = kozak_seq
            effect.kozak_strength = kozak_strength
            
            # Specific fields
            if 'aa_change' in effect_data:
                effect.aa_change = effect_data['aa_change']
            if 'alt_stop_exists' in effect_data:
                effect.alt_stop_exists = str(effect_data['alt_stop_exists'])
                effect.alt_stop_distance = str(effect_data.get('alt_stop_distance', 'NA'))
            if 'indel_length' in effect_data:
                ins_del = 'insertion' if effect_data.get('is_insertion') else 'deletion'
                effect.indel_type = f"{ins_del}_{effect_data['indel_length']}bp"
            
            effects_found.append(effect)
    
    return effects_found


# =============================================================================
# VCF Processing
# =============================================================================

def build_utr_regions_by_chrom(transcripts: Dict[str, TranscriptInfo],
                                chr_filter: Optional[Set[str]] = None) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Build a dictionary of 3'UTR regions organized by chromosome.
    
    This is used for efficient VCF region-based querying with tabix index.
    
    Args:
        transcripts: Dictionary of transcript information
        chr_filter: Optional set of chromosomes to include
        
    Returns:
        Dictionary mapping chromosome -> list of (start, end, transcript_id) tuples
    """
    utr_regions: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    
    for tid, tinfo in transcripts.items():
        # Apply chromosome filter
        if not chromosome_in_filter(tinfo.chrom, chr_filter):
            continue
        
        for utr_start, utr_end in tinfo.three_UTRs:
            utr_regions[tinfo.chrom].append((utr_start, utr_end, tid))
    
    # Sort regions by start position for each chromosome
    for chrom in utr_regions:
        utr_regions[chrom].sort(key=lambda x: x[0])
    
    return dict(utr_regions)


def merge_overlapping_regions(regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping genomic regions to reduce the number of tabix queries.
    
    Args:
        regions: List of (start, end) tuples
        
    Returns:
        List of merged non-overlapping (start, end) tuples
    """
    if not regions:
        return []
    
    # Sort by start position
    sorted_regions = sorted(regions, key=lambda x: x[0])
    merged = [sorted_regions[0]]
    
    for start, end in sorted_regions[1:]:
        last_start, last_end = merged[-1]
        # Merge if overlapping or adjacent (with 1bp gap tolerance)
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    
    return merged


def check_vcf_has_index(vcf_file: str) -> bool:
    """Check if VCF file has a tabix index (.tbi or .csi)."""
    tbi_path = vcf_file + '.tbi'
    csi_path = vcf_file + '.csi'
    return os.path.exists(tbi_path) or os.path.exists(csi_path)


def _process_region_batch(args: Tuple) -> List[List]:
    """
    Process a batch of genomic regions (for parallel execution).
    
    This function is designed to be called from ProcessPoolExecutor.
    Each process opens its own FASTA and VCF file to avoid sharing file handles.
    
    The batch can contain regions from one or multiple chromosomes, enabling
    parallel processing even when analyzing a single chromosome.
    
    Args:
        args: Tuple of (batch_id, region_list, vcf_file, fasta_path, transcripts_dict, 
                        dorfs_dict, csq_format, annotation_type, vcf_has_chr, utr_lookup)
              where region_list is [(chrom, start, end, utr_regions_for_region), ...]
    
    Returns:
        List of result rows (each row is a list for DataFrame)
    """
    (batch_id, region_list, vcf_file, fasta_path, transcripts, 
     dorfs, csq_format, annotation_type, vcf_has_chr, utr_lookup) = args
    
    results = []
    has_annotation = (annotation_type != 'none')
    
    try:
        # Each process opens its own FASTA and VCF
        fasta = pysam.FastaFile(fasta_path)
        vcf = VCF(vcf_file)
        
        for chrom, region_start, region_end, utr_regions_for_region in region_list:
            # Convert chromosome name to VCF format for querying
            vcf_chrom = normalize_chromosome(chrom, vcf_has_chr)
            region_str = f"{vcf_chrom}:{region_start}-{region_end}"
            try:
                for variant in vcf(region_str):
                    if variant.var_type not in ('snp', 'indel', 'mnp'):
                        continue
                    
                    pos = variant.POS
                    ref = variant.REF
                    
                    for alt_idx, alt in enumerate(variant.ALT):
                        if alt is None or alt == '*':
                            continue
                        alt = str(alt)
                        
                        symbol = ""
                        transcript_ids = []
                        
                        # Only try to extract annotation if VCF has annotation fields
                        if has_annotation:
                            # Extract annotation
                            symbol, transcript_ids = extract_annotation_from_variant(
                                variant, transcripts, csq_format
                            )
                        
                        # If no annotation found (or no annotation in VCF), check UTR overlap
                        if not transcript_ids:
                            for utr_start, utr_end, tid in utr_regions_for_region:
                                if utr_start <= pos <= utr_end:
                                    transcript_ids.append(tid)
                        
                        # Process each transcript
                        for tid in transcript_ids:
                            if tid not in transcripts:
                                continue
                            
                            transcript = transcripts[tid]
                            existing_dorfs = dorfs.get(tid, [])
                            gene_symbol = symbol if symbol else transcript.gene_name
                            
                            effects = process_variant(
                                chrom, pos, ref, alt, transcript, fasta, 
                                existing_dorfs, gene_symbol
                            )
                            
                            for effect in effects:
                                effect.transcript = tid
                                results.append(effect.to_list())
            except Exception as e:
                logger.warning(f"Error processing region {region_str}: {e}")
                continue
        
        fasta.close()
        vcf.close()
        
    except Exception as e:
        logger.error(f"Error in batch {batch_id}: {e}")
    
    return results


def _create_parallel_batches(utr_regions_by_chrom: Dict[str, List[Tuple[int, int, str]]], 
                              num_batches: int) -> List[List[Tuple]]:
    """
    Create balanced batches of regions for parallel processing.
    
    Distributes regions across batches to balance workload, regardless of
    chromosome boundaries. This enables parallel processing even for single
    chromosome analysis.
    
    Args:
        utr_regions_by_chrom: Dictionary of chrom -> [(start, end, tid), ...]
        num_batches: Target number of batches
        
    Returns:
        List of batches, each batch is [(chrom, start, end, relevant_utr_regions), ...]
    """
    # Collect all merged regions with their UTR info
    all_regions = []
    
    for chrom, regions_with_tids in utr_regions_by_chrom.items():
        # Get just (start, end) for merging
        regions = [(r[0], r[1]) for r in regions_with_tids]
        merged = merge_overlapping_regions(regions)
        
        for start, end in merged:
            # Find UTRs that overlap this merged region
            relevant_utrs = [(s, e, t) for s, e, t in regions_with_tids 
                            if not (e < start or s > end)]
            all_regions.append((chrom, start, end, relevant_utrs))
    
    if not all_regions:
        return []
    
    # Sort by region size (descending) for better load balancing
    all_regions.sort(key=lambda x: x[2] - x[1], reverse=True)
    
    # Distribute regions to batches using round-robin
    batches = [[] for _ in range(min(num_batches, len(all_regions)))]
    for i, region in enumerate(all_regions):
        batches[i % len(batches)].append(region)
    
    return batches


def process_vcf_file(vcf_file: str, transcripts: Dict[str, TranscriptInfo],
                    dorfs: Dict[str, List[dORFEntry]], fasta: pysam.FastaFile,
                    output_prefix: str, threads: int = 1,
                    chr_filter: Optional[Set[str]] = None,
                    region_filter: Optional[GenomicRegion] = None,
                    fasta_path: Optional[str] = None) -> int:
    """
    Process VCF file and annotate variants.
    
    Performance optimizations:
    1. When VCF is indexed (tabix), uses region-based querying to only fetch 
       variants in 3'UTR regions.
    2. When threads > 1 and VCF is indexed, uses parallel processing by chromosome.
    
    Args:
        vcf_file: Path to VCF file
        transcripts: Dictionary of transcript information
        dorfs: Dictionary of known dORFs by transcript
        fasta: Open pysam FastaFile (used for single-threaded mode)
        output_prefix: Output file prefix
        threads: Number of threads for parallel processing
        chr_filter: Optional set of chromosomes to process
        region_filter: Optional genomic region to filter
        fasta_path: Path to FASTA file (required for parallel mode)
    
    Returns:
        Number of effects found
    """
    
    filter_msg = ""
    if chr_filter or region_filter:
        filter_msg = " (with chromosome/region filter)"
    logger.info(f"Processing VCF file: {vcf_file}{filter_msg}")
    
    all_results = []
    processed = 0
    skipped_by_filter = 0
    
    # Build position-to-transcript index for faster 3'UTR overlap lookup
    utr_regions_by_chrom = build_utr_regions_by_chrom(transcripts, chr_filter)
    
    # Check if VCF has index for region-based querying
    has_index = check_vcf_has_index(vcf_file)
    use_region_query = has_index and not region_filter  # Don't use with explicit region filter
    
    # Determine if parallel processing is possible
    # Parallel mode works even for single chromosome by splitting regions into batches
    use_parallel = (threads > 1 and has_index and fasta_path and 
                    len(utr_regions_by_chrom) >= 1 and not region_filter)
    
    if use_region_query:
        logger.info("VCF index detected - using region-based querying for 3'UTR regions")
    if use_parallel:
        logger.info(f"Parallel processing enabled with {threads} workers")
    
    vcf = VCF(vcf_file)
    
    # Detect CSQ format from VCF header for VEP parsing
    vcf_header = vcf.raw_header if hasattr(vcf, 'raw_header') else ""
    csq_format = detect_csq_format(vcf_header)
    
    # Detect annotation type
    annotation_type = detect_vcf_annotation_type(vcf_header)
    has_annotation = (annotation_type != 'none')
    
    # Detect VCF chromosome naming format
    vcf_has_chr = detect_vcf_chr_format(vcf_file)
    logger.info(f"VCF chromosome format: {'chr prefix' if vcf_has_chr else 'no chr prefix'}")
    
    # Log annotation detection results
    if annotation_type == 'custom_vep':
        logger.info("Detected custom VEP annotation fields in VCF")
    elif annotation_type == 'standard_vep':
        logger.info(f"Detected standard VEP CSQ annotation: {csq_format[:50] if csq_format else 'format unknown'}...")
    elif annotation_type == 'annovar':
        logger.info("Detected ANNOVAR annotation fields in VCF")
    else:
        logger.info("No VEP/ANNOVAR/custom annotation detected - using position-based annotation")
        logger.info("Will annotate variants based on position, GTF, FASTA and dORF files directly")
    
    def process_single_variant(variant, chrom: str, pos: int, ref: str) -> None:
        """Process a single variant and append results to all_results."""
        nonlocal processed
        
        # Handle multi-allelic sites
        for alt_idx, alt in enumerate(variant.ALT):
            if alt is None or alt == '*':
                continue
            
            alt = str(alt)
            symbol = ""
            transcript_ids = []
            
            # Only try to extract annotation if VCF has annotation fields
            if has_annotation:
                # Extract annotation using unified function (supports custom VEP, standard VEP CSQ, ANNOVAR)
                symbol, transcript_ids = extract_annotation_from_variant(variant, transcripts, csq_format)
            
            # If no annotation found (or no annotation in VCF), check UTR regions for this chromosome
            if not transcript_ids:
                # Use pre-built index for faster lookup
                if chrom in utr_regions_by_chrom:
                    for utr_start, utr_end, tid in utr_regions_by_chrom[chrom]:
                        if utr_start <= pos <= utr_end:
                            transcript_ids.append(tid)
                else:
                    # Try with normalized chromosome name
                    norm_chrom = normalize_chromosome(chrom, False)
                    alt_chrom = normalize_chromosome(chrom, True)
                    check_chroms = [norm_chrom, alt_chrom]
                    
                    for check_chrom in check_chroms:
                        if check_chrom in utr_regions_by_chrom:
                            for utr_start, utr_end, tid in utr_regions_by_chrom[check_chrom]:
                                if utr_start <= pos <= utr_end:
                                    transcript_ids.append(tid)
                            break
            
            # Process each affected transcript
            for tid in transcript_ids:
                if tid not in transcripts:
                    continue
                
                transcript = transcripts[tid]
                existing_dorfs = dorfs.get(tid, [])
                gene_symbol = symbol if symbol else transcript.gene_name
                
                effects = process_variant(
                    chrom, pos, ref, alt, transcript, fasta, existing_dorfs, gene_symbol
                )
                
                for effect in effects:
                    effect.transcript = tid
                    all_results.append(effect.to_list())
        
        processed += 1
        if processed % 10000 == 0:
            logger.info(f"Processed {processed} variants, found {len(all_results)} effects...")
    
    if use_parallel:
        # =====================================================================
        # Parallel processing mode: process region batches in parallel
        # This works even for single chromosome by splitting regions into batches
        # =====================================================================
        vcf.close()  # Close VCF, each worker will open its own
        
        # Create balanced batches of regions
        batches = _create_parallel_batches(utr_regions_by_chrom, threads * 2)
        
        if not batches:
            logger.warning("No regions to process in parallel mode")
        else:
            # Count total regions for logging
            total_regions = sum(len(batch) for batch in batches)
            logger.info(f"Splitting {total_regions} regions into {len(batches)} batches "
                       f"for {threads} parallel workers")
            
            # Build UTR lookup for fast access in workers
            utr_lookup = utr_regions_by_chrom
            
            # Prepare tasks
            tasks = []
            for batch_id, batch in enumerate(batches):
                tasks.append((
                    batch_id,
                    batch,
                    vcf_file,
                    fasta_path,
                    transcripts,
                    dorfs,
                    csq_format,
                    annotation_type,
                    vcf_has_chr,
                    utr_lookup
                ))
            
            # Execute in parallel
            with ProcessPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(_process_region_batch, task): task[0] 
                          for task in tasks}
                
                completed = 0
                for future in as_completed(futures):
                    batch_id = futures[future]
                    try:
                        batch_results = future.result()
                        all_results.extend(batch_results)
                        completed += 1
                        logger.info(f"Completed batch {batch_id + 1}/{len(batches)}: "
                                   f"{len(batch_results)} effects")
                    except Exception as e:
                        logger.error(f"Error processing batch {batch_id}: {e}")
        
        processed = len(all_results)  # Approximate
        
    elif use_region_query:
        # =====================================================================
        # Single-threaded region-based querying
        # =====================================================================
        for chrom, regions_with_tids in utr_regions_by_chrom.items():
            # Apply chromosome filter (already filtered in build_utr_regions_by_chrom, but double-check)
            if not chromosome_in_filter(chrom, chr_filter):
                continue
            
            # Extract just (start, end) and merge overlapping regions
            regions = [(r[0], r[1]) for r in regions_with_tids]
            merged_regions = merge_overlapping_regions(regions)
            
            logger.debug(f"Querying {len(merged_regions)} merged 3'UTR regions on {chrom}")
            
            # Convert chromosome name to VCF format for querying
            vcf_chrom = normalize_chromosome(chrom, vcf_has_chr)
            
            for region_start, region_end in merged_regions:
                # cyvcf2 region query format: "chrom:start-end"
                # Use VCF's chromosome naming format
                region_str = f"{vcf_chrom}:{region_start}-{region_end}"
                try:
                    for variant in vcf(region_str):
                        # Skip non-SNV/indel
                        if variant.var_type not in ('snp', 'indel', 'mnp'):
                            continue
                        
                        # Convert VCF chromosome name to FASTA format (our internal standard)
                        variant_chrom = standardize_chromosome_name(variant.CHROM, config.chr_prefix)
                        process_single_variant(variant, variant_chrom, variant.POS, variant.REF)
                except Exception as e:
                    logger.warning(f"Region query failed for {region_str}: {e}")
                    continue
        
        vcf.close()
    else:
        # =====================================================================
        # Traditional full-file iteration (no index)
        # =====================================================================
        for variant in vcf:
            # Skip non-SNV/indel
            if variant.var_type not in ('snp', 'indel', 'mnp'):
                continue
            
            # Convert VCF chromosome name to FASTA format (our internal standard)
            chrom = standardize_chromosome_name(variant.CHROM, config.chr_prefix)
            pos = variant.POS
            ref = variant.REF
            
            # Apply chromosome filter
            if not chromosome_in_filter(chrom, chr_filter):
                skipped_by_filter += 1
                continue
            
            # Apply region filter
            if region_filter and not region_filter.contains(chrom, pos):
                skipped_by_filter += 1
                continue
            
            process_single_variant(variant, chrom, pos, ref)
        
        vcf.close()
    
    # Log cache statistics (only meaningful for single-threaded mode)
    if not use_parallel:
        cache_stats = _sequence_cache.stats()
        logger.info(f"Sequence cache stats: {cache_stats['hits']} hits, {cache_stats['misses']} misses "
                    f"({cache_stats['hit_rate']} hit rate)")
    
    logger.info(f"Total: processed {processed} variants, found {len(all_results)} effects")
    
    # Write results
    if all_results:
        df = pd.DataFrame(all_results, columns=OUTPUT_COLUMNS)
        output_file = f"{output_prefix}_dORFannotator.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Results written to: {output_file}")
    else:
        logger.warning("No dORF effects found")
    
    return len(all_results)


# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main entry point."""
    
    parser = argparse.ArgumentParser(
        description='dORFannotator - Annotate 3\'UTR variants affecting downstream ORFs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full analysis
  python dORFannotator.py --genome hg38.fa --vcf input.vcf.gz --dorf known_dORF.bed --annotation gencode.gff3.gz --out results

  # Filter by single chromosome (faster for testing)
  python dORFannotator.py --genome hg38.fa --vcf input.vcf.gz --dorf known_dORF.bed --annotation gencode.gff3.gz --out results --chr chr1

  # Filter by multiple chromosomes
  python dORFannotator.py --genome hg38.fa --vcf input.vcf.gz --dorf known_dORF.bed --annotation gencode.gff3.gz --out results --chr chr1,chr2,chrX

  # Filter by specific region
  python dORFannotator.py --genome hg38.fa --vcf input.vcf.gz --dorf known_dORF.bed --annotation gencode.gff3.gz --out results --region chr1:1000000-2000000

Report bugs to: liusihan@wchscu.cn
Homepage: https://github.com/liusihan/dORFannotator
        """
    )
    
    # Required arguments
    parser.add_argument('--genome', '-g', required=True,
                       help='Reference genome FASTA file (with .fai index)')
    parser.add_argument('--vcf', '-v', required=True,
                       help='Input VCF file (bgzipped with .tbi/.csi index)')
    parser.add_argument('--dorf', '-d', required=True,
                       help='Known dORF BED file')
    parser.add_argument('--annotation', '-a', required=True,
                       help='Gene annotation GTF/GFF file')
    parser.add_argument('--out', '-o', required=True,
                       help='Output file prefix')
    
    # Optional arguments - Filtering
    parser.add_argument('--chr', type=str, default=None,
                       help='Filter by chromosome (e.g., chr1 or 1). Multiple chromosomes separated by comma (e.g., chr1,chr2,chr3)')
    parser.add_argument('--region', type=str, default=None,
                       help='Filter by genomic region (e.g., chr1:1000000-2000000)')
    
    # Optional arguments - Processing
    parser.add_argument('--threads', '-t', type=int, default=1,
                       help='Number of threads (default: 1)')
    parser.add_argument('--start-codons', type=str, default=None,
                       help='Comma-separated start codons (default: ATG,CTG,GTG,TTG,ACG)')
    parser.add_argument('--stop-codons', type=str, default=None,
                       help='Comma-separated stop codons (default: TAG,TAA,TGA)')
    parser.add_argument('--min-length', type=int, default=30,
                       help='Minimum ORF length in bp (default: 30)')
    parser.add_argument('--max-length', type=int, default=303,
                       help='Maximum ORF length in bp (default: 303)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--version', action='version', version=f'dORFannotator {__version__}')
    
    args = parser.parse_args()
    
    # Setup logging
    global logger
    logger = setup_logging(args.verbose)
    
    # Print header
    print("\n" + "="*70)
    print(f"  dORFannotator v{__version__}")
    print("  Annotate 3'UTR variants affecting downstream ORFs")
    print("  (C) 2024 Sihan Liu - West China Hospital")
    print("="*70 + "\n")
    
    start_time = time.time()
    
    # Validate input files
    if not validate_input_files(args.genome, args.vcf, args.annotation, args.dorf):
        sys.exit(1)
    
    # Configure parameters
    if args.start_codons:
        config.start_codons = set(args.start_codons.upper().split(','))
    if args.stop_codons:
        config.stop_codons = set(args.stop_codons.upper().split(','))
    config.min_orf_length = args.min_length
    config.max_orf_length = args.max_length
    config.threads = args.threads
    
    logger.info(f"Start codons: {config.start_codons}")
    logger.info(f"Stop codons: {config.stop_codons}")
    logger.info(f"ORF length range: {config.min_orf_length}-{config.max_orf_length} bp")
    
    # Parse chromosome and region filters
    chr_filter = parse_chromosome_filter(args.chr)
    region_filter = parse_region_string(args.region) if args.region else None
    
    if chr_filter:
        logger.info(f"Chromosome filter: {', '.join(sorted(list(chr_filter)[:10]))}{'...' if len(chr_filter) > 10 else ''}")
    if region_filter:
        region_str = f"{region_filter.chrom}"
        if region_filter.start or region_filter.end:
            region_str += f":{region_filter.start or ''}-{region_filter.end or ''}"
        logger.info(f"Region filter: {region_str}")
    
    # Load input files
    try:
        # Check chromosome format consistency across all files
        chr_formats = check_chromosome_format_consistency(
            args.genome, args.vcf, args.annotation, args.dorf
        )
        
        # Use FASTA format as our internal standard
        config.chr_prefix = chr_formats['fasta']
        logger.info(f"\nUsing FASTA chromosome format as internal standard: "
                   f"{'chr prefix' if config.chr_prefix else 'no chr prefix'}")
        logger.info(f"All files will be standardized to this format during processing.\n")
        
        # Parse GTF with chromosome standardization
        transcripts = parse_gtf_file(args.annotation, chr_filter, region_filter, 
                                     target_chr_format=config.chr_prefix)
        
        # Parse dORF BED with chromosome standardization
        dorfs = parse_dorf_bed(args.dorf, chr_filter, region_filter,
                              target_chr_format=config.chr_prefix)
        
        # Open genome FASTA
        fasta = pysam.FastaFile(args.genome)
        
        # Process VCF
        # Pass fasta_path for parallel processing support
        effect_count = process_vcf_file(
            args.vcf, transcripts, dorfs, fasta, args.out, args.threads,
            chr_filter, region_filter, fasta_path=args.genome
        )
        
        fasta.close()
        
        elapsed = time.time() - start_time
        logger.info(f"Analysis completed in {elapsed:.2f} seconds")
        logger.info(f"Found {effect_count} dORF effects")
        
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print("\n" + "="*70)
    print("  Analysis complete!")
    print("  Report bugs to: liusihan@wchscu.cn")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()


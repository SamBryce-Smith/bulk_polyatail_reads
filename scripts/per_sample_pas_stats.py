#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import NamedTuple
import sys

# see https://stackoverflow.com/a/50038614 - suggested syntax for python 3.6+ when want type hints
class Sample(NamedTuple):
    name: str
    stats_path: Path
    bed_path: Path

class Results(NamedTuple):
    sample_name: str
    patr_count: int
    primary_count: int
    patr_percentage: float
    patr_per_million: float

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Analyze poly-A tail read alignments.')
    parser.add_argument('input_tsv', type=Path, help='Input TSV file with sample information')
    parser.add_argument('output_tsv', type=Path, help='Output TSV file for results')
    parser.add_argument('--percentage_dp', type=int, default=6,
                       help='Number of decimal places for percentage values (default: 6)')
    return parser.parse_args()

def validate_paths(sample: Sample):
    """Check if the paths for a sample exist."""
    if not sample.stats_path.exists():
        raise FileNotFoundError(
            f"Stats file not found for sample '{sample.name}': {sample.stats_path}"
        )
    if not sample.bed_path.exists():
        raise FileNotFoundError(
            f"BED file not found for sample '{sample.name}': {sample.bed_path}"
        )

def get_primary_alignments(stats_path: Path) -> int:
    """Extract primary alignment count from stats TSV."""
    with open(stats_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'stat' not in reader.fieldnames or 'value' not in reader.fieldnames:
            raise ValueError(f"Stats file missing required columns 'stat' and 'value': {stats_path}")

        for row in reader:
            if row['stat'] == 'primary_alignments':
                return int(row['value'])
    raise ValueError(f"No primary_alignments found in {stats_path}")

def get_patr_count(bed_path: Path) -> int:
    """Sum the score column (field 5) from BED6 file."""
    total = 0
    with open(bed_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            total += int(fields[4])
    return total

def analyze_sample(sample: Sample) -> Results:
    """Analyze a single sample and return results."""
    primary_count = get_primary_alignments(sample.stats_path)
    patr_count = get_patr_count(sample.bed_path)
    patr_percentage = (patr_count / primary_count * 100) if primary_count > 0 else 0.0
    # Calculate reads per million (RPM)
    # 1. Calculate the "per million" scaling factor
    scaling_factor = primary_count / 1_000_000 if primary_count > 0 else 0.0
    # 2. Normalize counts by scaling factor to get RPM
    patr_per_million = patr_count / scaling_factor if scaling_factor > 0 else 0.0
    
    
    return Results(
        sample_name=sample.name,
        patr_count=patr_count,
        primary_count=primary_count,
        patr_percentage=patr_percentage,
        patr_per_million=patr_per_million
    )

def process_samples(input_tsv: Path, output_tsv: Path, percentage_dp: int):
    """Process samples and write results directly to output file."""
    # Create output file and write header
    with open(output_tsv, 'w') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(['sample_name', 'patr_count', 'primary_count', 'patr_percentage', 'patr_per_million'])
        
        # Process samples one at a time
        with open(input_tsv) as in_f:
            reader = csv.DictReader(in_f, delimiter='\t')
            
            # Verify we have the required columns
            required_columns = {'sample_name', 'stats_path', 'bed_path'}
            missing_columns = required_columns - set(reader.fieldnames)
            if missing_columns:
                raise ValueError(f"Input TSV missing required columns: {missing_columns}")
            
            for row in reader:
                sample = Sample(
                    name=row['sample_name'],
                    stats_path=Path(row['stats_path']),
                    bed_path=Path(row['bed_path'])
                )
                
                # Validate paths before processing
                validate_paths(sample)
                
                # Process sample and write results immediately
                result = analyze_sample(sample)
                writer.writerow([
                    result.sample_name,
                    result.patr_count,
                    result.primary_count,
                    f"{result.patr_percentage:.{percentage_dp}f}",
                    f"{result.patr_per_million:.{percentage_dp}f}"
                ])
                        

def main():
    args = parse_args()
    
    try:
        process_samples(args.input_tsv, args.output_tsv, args.percentage_dp)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
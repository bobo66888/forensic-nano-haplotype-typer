#!/usr/bin/env python3
"""
SAM File SNP Analysis Tool with ADR Calculation

This script analyzes SAM files from third-generation sequencing data,
extracts SNP patterns from aligned reads, calculates Allelic Depth Ratio (ADR),
and determines genotype (heterozygous/homozygous) based on ADR threshold.

Author: Bioinformatics Tool
Version: 2.0
"""

import os
import sys
import pandas as pd
from collections import Counter
import gzip
import re


def parse_cigar(cigar_string):
    """
    Parse CIGAR string and return list of operations.

    Args:
        cigar_string: CIGAR string from SAM file

    Returns:
        list: List of tuples (length, operation)
    """
    operations = []
    current_num = ""

    for char in cigar_string:
        if char.isdigit():
            current_num += char
        else:
            if current_num:
                operations.append((int(current_num), char))
                current_num = ""

    return operations


def find_base_at_position(read_seq, cigar_ops, read_start, target_pos):
    """
    Find base at specific reference position considering CIGAR operations.

    Args:
        read_seq: Read sequence
        cigar_ops: Parsed CIGAR operations
        read_start: Starting position of read in reference
        target_pos: Target reference position

    Returns:
        str: Base at target position or 'N'/'D' for unavailable/deletion
    """
    ref_pos = read_start
    read_idx = 0

    for op_length, op_type in cigar_ops:
        if op_type in "M=X":
            op_end = ref_pos + op_length - 1
            if ref_pos <= target_pos <= op_end:
                offset = target_pos - ref_pos
                if read_idx + offset < len(read_seq):
                    return read_seq[read_idx + offset]
                else:
                    return 'N'
            ref_pos += op_length
            read_idx += op_length

        elif op_type == "I":
            read_idx += op_length

        elif op_type in "DN":
            op_end = ref_pos + op_length - 1
            if ref_pos <= target_pos <= op_end:
                return 'D'
            ref_pos += op_length

        elif op_type == "S":
            read_idx += op_length

        elif op_type == "H":
            continue

        elif op_type == "P":
            continue

        else:
            continue

    return 'N'


def parse_sam_line(line):
    """
    Parse a single line from SAM file.

    Args:
        line: Line from SAM file

    Returns:
        dict: Dictionary with SAM fields or None if malformed
    """
    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None

    return {
        'qname': fields[0],
        'flag': int(fields[1]),
        'rname': fields[2],
        'pos': int(fields[3]),
        'mapq': int(fields[4]),
        'cigar': fields[5],
        'rnext': fields[6],
        'pnext': int(fields[7]),
        'tlen': int(fields[8]),
        'seq': fields[9],
        'qual': fields[10]
    }


def get_bases_at_positions(read_data, positions):
    """
    Extract bases at multiple SNP positions from a read.

    Args:
        read_data: Dictionary with read information
        positions: List of SNP positions in reference

    Returns:
        str: Concatenated bases at SNP positions
    """
    cigar_ops = parse_cigar(read_data['cigar'])
    bases = []

    for pos in positions:
        base = find_base_at_position(
            read_data['seq'],
            cigar_ops,
            read_data['pos'],
            pos
        )
        bases.append(base)

    return ''.join(bases)


def read_snp_positions_from_excel(excel_file):
    """
    Read SNP position information from Excel file.

    Args:
        excel_file: Path to Excel file with SNP positions

    Returns:
        dict: Dictionary with site names as keys and SNP positions as values
    """
    try:
        # Read Excel file (assuming Chinese column names from original data)
        df = pd.read_excel(excel_file, sheet_name=0)

        # Dictionary to store SNP positions
        snp_positions = {}

        # Iterate through each row
        for _, row in df.iterrows():
            # Handle site name (assuming Chinese column name 'Locus')
            site_name = str(row['Locus']).strip()
            snp_string = row['SNP_position']

            # Process SNP position string
            if pd.isna(snp_string):
                positions = []
            else:
                # Convert string to list of integers
                positions = []
                snp_string = str(snp_string)
                # Support comma, space, semicolon separators
                for item in re.split(r'[,;\s]+', snp_string):
                    if item.strip().isdigit():
                        positions.append(int(item.strip()))

            snp_positions[site_name] = positions

        print(f"Successfully read SNP information for {len(snp_positions)} sites from Excel file")

        return snp_positions

    except Exception as e:
        print(f"Error reading Excel file: {e}")
        import traceback
        traceback.print_exc()
        return {}


def get_target_reference_from_sam(sam_file):
    """
    Extract the actual aligned target reference from SAM file.

    Args:
        sam_file: Path to SAM file

    Returns:
        tuple: (target_reference_name, reference_length) or (None, None)
    """
    # Determine if file is compressed
    is_gzipped = sam_file.endswith('.gz')

    # Open file
    if is_gzipped:
        file_opener = gzip.open
        mode = 'rt'
    else:
        file_opener = open
        mode = 'r'

    print("Identifying target reference sequence from SAM file...")

    target_ref = None
    ref_length = None

    with file_opener(sam_file, mode) as f:
        # First scan header to get all reference information
        references = {}
        for line in f:
            if line.startswith('@SQ'):
                fields = line.strip().split('\t')
                ref_name = None
                ref_len = None

                for field in fields:
                    if field.startswith('SN:'):
                        ref_name = field[3:]  # Remove 'SN:'
                    elif field.startswith('LN:'):
                        ref_len = int(field[3:])  # Remove 'LN:'

                if ref_name and ref_len is not None:
                    references[ref_name] = ref_len

            # Find first alignment line
            elif not line.startswith('@'):
                # Parse alignment line
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    rname = fields[2]
                    if rname != '*':  # Ensure aligned read
                        target_ref = rname
                        break

        # If target found, get its length
        if target_ref and target_ref in references:
            ref_length = references[target_ref]

    if target_ref:
        print(f"Target reference identified: {target_ref} (length: {ref_length})")
    else:
        print("Warning: Could not identify target reference sequence from SAM file")

    return target_ref, ref_length


def find_matching_site_in_excel(target_ref, excel_snp_positions):
    """
    Find matching site in Excel data based on target reference.

    Args:
        target_ref: Target reference name from SAM file
        excel_snp_positions: Dictionary of SNP positions from Excel

    Returns:
        tuple: (matched_excel_site_name, SNP_positions_list) or (None, None)
    """
    print(f"Looking for site matching reference '{target_ref}' in Excel data...")

    # Strategy 1: Direct match
    if target_ref in excel_snp_positions:
        print(f"  Direct match successful: '{target_ref}'")
        return target_ref, excel_snp_positions[target_ref]

    # Strategy 2: Match after removing 'S' prefix
    if target_ref.startswith('S'):
        target_without_s = target_ref[1:]
        if target_without_s in excel_snp_positions:
            print(f"  Match after removing 'S' prefix: '{target_ref}' -> '{target_without_s}'")
            return target_without_s, excel_snp_positions[target_without_s]

    # Strategy 3: Partial match (extract number part)
    match = re.search(r'S?(\d+-\w+)', target_ref)
    if match:
        possible_name = match.group(1)
        if possible_name in excel_snp_positions:
            print(f"  Partial match successful: '{target_ref}' -> '{possible_name}'")
            return possible_name, excel_snp_positions[possible_name]

    # Strategy 4: Match by identifier (e.g., MH135189)
    parts = target_ref.split('-')
    if len(parts) > 1:
        # Take last part as possible identifier
        identifier = parts[-1]
        for excel_site, positions in excel_snp_positions.items():
            if identifier in excel_site:
                print(f"  Identifier match successful: '{target_ref}' -> '{excel_site}'")
                return excel_site, positions

    print(f"  No matching site found")
    return None, None


def calculate_adr_and_genotype(pattern_counts, adr_threshold=0.1):
    """
    Calculate Allelic Depth Ratio (ADR) and determine genotype.

    ADR = count of second most frequent haplotype / count of most frequent haplotype

    Args:
        pattern_counts: Counter object with haplotypes and their counts
        adr_threshold: ADR threshold for genotype determination (default: 0.1)

    Returns:
        dict: Dictionary with ADR analysis results
    """
    results = {
        'genotype': 'unknown',
        'adr_value': 0,
        'haplotype1': None,
        'count1': 0,
        'haplotype2': None,
        'count2': 0,
        'note': ''
    }

    # Check if there is data
    if not pattern_counts or len(pattern_counts) == 0:
        results['note'] = 'No covered reads'
        return results

    # Get most common haplotypes
    most_common = pattern_counts.most_common()

    if len(most_common) == 0:
        results['note'] = 'No valid haplotypes'
        return results

    # Get top haplotype
    haplotype1, count1 = most_common[0]
    results['haplotype1'] = haplotype1
    results['count1'] = count1

    # If there is a second haplotype
    if len(most_common) > 1:
        haplotype2, count2 = most_common[1]
        results['haplotype2'] = haplotype2
        results['count2'] = count2

        # Calculate ADR
        if count1 > 0:
            adr = count2 / count1
            results['adr_value'] = adr

            # Determine genotype
            if adr_threshold <= adr <= 1:
                results['genotype'] = 'heterozygous'
                results['note'] = f'Heterozygous (ADR={adr:.3f})'
            elif 0 <= adr < adr_threshold:
                results['genotype'] = 'homozygous'
                # For homozygous, treat top haplotype as two identical haplotypes
                results['haplotype2'] = haplotype1
                results['count2'] = count1 / 2
                results['count1'] = count1 / 2
                results['note'] = f'Homozygous (ADR={adr:.3f})'
            else:
                results['genotype'] = 'unknown'
                results['note'] = f'Abnormal ADR value (ADR={adr:.3f})'
    else:
        # Only one haplotype detected, treat as homozygous
        results['genotype'] = 'homozygous'
        results['haplotype2'] = haplotype1
        results['count2'] = count1 / 2
        results['count1'] = count1 / 2
        results['note'] = 'Only one haplotype detected, treated as homozygous'

    return results


def analyze_single_sam_file(sam_file, excel_file, output_file, adr_output_file=None, min_mapq=20, adr_threshold=0.1):
    """
    Analyze single SAM file for target reference sequence.

    Args:
        sam_file: Path to SAM file
        excel_file: Path to Excel file with SNP positions
        output_file: Path for main output Excel file
        adr_output_file: Path for ADR analysis output Excel file (optional)
        min_mapq: Minimum mapping quality threshold (default: 20)
        adr_threshold: ADR threshold for genotype determination (default: 0.1)
    """
    # Check if files exist
    if not os.path.exists(sam_file):
        print(f"Error: SAM file '{sam_file}' does not exist")
        return

    if not os.path.exists(excel_file):
        print(f"Error: Excel file '{excel_file}' does not exist")
        return

    print("=" * 60)
    print("Starting SAM file analysis")
    print("=" * 60)

    # Step 1: Read SNP positions from Excel
    print("\n1. Reading SNP position information from Excel file...")
    excel_snp_positions = read_snp_positions_from_excel(excel_file)

    if not excel_snp_positions:
        print("Error: Could not read SNP position information from Excel file")
        return

    # Step 2: Determine target reference from SAM file
    print("\n2. Determining target reference sequence from SAM file...")
    target_ref, ref_length = get_target_reference_from_sam(sam_file)

    if not target_ref:
        print("Error: Could not determine target reference sequence from SAM file")
        return

    # Step 3: Find matching site in Excel
    print("\n3. Finding matching site in Excel...")
    excel_site, snp_positions = find_matching_site_in_excel(target_ref, excel_snp_positions)

    if not excel_site or not snp_positions:
        print(f"Error: Reference '{target_ref}' has no corresponding SNP information in Excel")
        print(f"Excel sites: {list(excel_snp_positions.keys())[:10]}")  # Show first 10
        if len(excel_snp_positions) > 10:
            print(f"  ... Total {len(excel_snp_positions)} sites")
        return

    print(f"\nMatch successful:")
    print(f"  SAM reference: {target_ref}")
    print(f"  Excel site: {excel_site}")
    print(f"  SNP positions: {snp_positions}")
    print(f"  Number of SNPs: {len(snp_positions)}")

    # Step 4: Analyze reads in SAM file
    print("\n4. Analyzing reads in SAM file...")

    # Initialize statistics
    pattern_counts = Counter()
    total_reads = 0
    covering_reads = 0
    skipped_low_quality = 0
    unmapped_reads = 0
    other_ref_reads = 0

    # Determine if file is compressed
    is_gzipped = sam_file.endswith('.gz')

    # Open SAM file
    if is_gzipped:
        file_opener = gzip.open
        mode = 'rt'
    else:
        file_opener = open
        mode = 'r'

    read_count = 0
    processed_count = 0

    with file_opener(sam_file, mode) as f:
        for line_num, line in enumerate(f, 1):
            # Skip header lines
            if line.startswith('@'):
                continue

            read_count += 1

            # Show progress every 10000 lines
            if read_count % 10000 == 0:
                print(f"  Processed {read_count} reads...")

            # Parse SAM line
            read_data = parse_sam_line(line)
            if not read_data:
                continue

            # Get reference name
            rname = read_data['rname']

            # Count unmapped reads
            if rname == '*':
                unmapped_reads += 1
                continue

            # Only process reads aligned to target reference
            if rname != target_ref:
                other_ref_reads += 1
                continue

            total_reads += 1
            processed_count += 1

            # Skip low quality reads
            if read_data['mapq'] < min_mapq:
                skipped_low_quality += 1
                continue

            # Get bases at all target positions
            base_pattern = get_bases_at_positions(read_data, snp_positions)

            # Check if all positions have valid bases (not 'N')
            if 'N' not in base_pattern:
                covering_reads += 1
                pattern_counts[base_pattern] += 1

    print(f"\nRead processing completed:")
    print(f"  Total reads: {read_count}")
    print(f"  Reads aligned to target reference: {total_reads}")
    print(f"  Successfully processed reads: {processed_count}")
    print(f"  Reads aligned to other references: {other_ref_reads}")
    print(f"  Unmapped reads: {unmapped_reads}")

    # Step 5: Calculate ADR and determine genotype
    print(f"\n5. Calculating ADR and determining genotype (threshold: {adr_threshold})...")
    adr_results = calculate_adr_and_genotype(pattern_counts, adr_threshold)

    # Step 6: Generate main output results
    print("\n6. Generating main output results...")

    # Create output DataFrame
    all_results = []

    # Add basic information
    all_results.append({
        'Site': excel_site,
        'Information': 'SNP positions',
        'Value': str(snp_positions)
    })

    # Add statistics
    if total_reads > 0:
        coverage_rate = covering_reads / total_reads * 100
    else:
        coverage_rate = 0

    all_results.append({
        'Site': excel_site,
        'Information': 'Total reads',
        'Value': total_reads
    })

    all_results.append({
        'Site': excel_site,
        'Information': 'Covering reads',
        'Value': covering_reads
    })

    all_results.append({
        'Site': excel_site,
        'Information': 'Coverage rate',
        'Value': f"{coverage_rate:.2f}%"
    })

    all_results.append({
        'Site': excel_site,
        'Information': 'Low quality filtered',
        'Value': skipped_low_quality
    })

    all_results.append({
        'Site': excel_site,
        'Information': 'ADR analysis',
        'Value': adr_results['note']
    })

    all_results.append({
        'Site': excel_site,
        'Information': 'ADR value',
        'Value': f"{adr_results['adr_value']:.4f}"
    })

    # Add separator
    all_results.append({
        'Site': '',
        'Information': '',
        'Value': ''
    })

    # Add haplotype statistics header
    all_results.append({
        'Site': excel_site,
        'Information': 'Haplotype pattern',
        'Value': 'Count'
    })

    # Add each haplotype pattern statistics
    if pattern_counts:
        for pattern, count in pattern_counts.most_common():
            percentage = count / covering_reads * 100 if covering_reads > 0 else 0
            all_results.append({
                'Site': excel_site,
                'Information': pattern,
                'Value': f"{count} ({percentage:.2f}%)"
            })
    else:
        all_results.append({
            'Site': excel_site,
            'Information': 'No covering reads',
            'Value': 0
        })

    # Create DataFrame
    df_output = pd.DataFrame(all_results)

    # Save to Excel file
    try:
        df_output.to_excel(output_file, index=False)
        print(f"\nMain results saved to: {output_file}")

    except Exception as e:
        print(f"Error saving main results to Excel file: {e}")
        import traceback
        traceback.print_exc()

    # Step 7: Generate ADR analysis output file
    if adr_output_file:
        print("\n7. Generating ADR analysis results file...")

        # Create ADR results DataFrame
        adr_data = {
            'Site': [excel_site],
            'Genotype': [adr_results['genotype']],
            'ADR': [f"{adr_results['adr_value']:.4f}"],
            'Haplotype1': [adr_results['haplotype1']],
            'Count1': [adr_results['count1']],
            'Haplotype2': [adr_results['haplotype2']],
            'Count2': [adr_results['count2']]
        }

        # Create DataFrame
        df_adr = pd.DataFrame(adr_data)

        # Save to Excel file
        try:
            df_adr.to_excel(adr_output_file, index=False)
            print(f"ADR analysis results saved to: {adr_output_file}")

        except Exception as e:
            print(f"Error saving ADR analysis results to Excel file: {e}")
            import traceback
            traceback.print_exc()

    # Step 8: Print summary
    print("\n" + "=" * 60)
    print("Analysis Summary:")
    print("=" * 60)
    print(f"Target reference: {target_ref}")
    print(f"Corresponding Excel site: {excel_site}")
    print(f"SNP positions: {snp_positions}")
    print(f"Number of SNPs: {len(snp_positions)}")
    print(f"Total reads: {total_reads}")
    print(f"Covering reads: {covering_reads}")
    print(f"Coverage rate: {coverage_rate:.2f}%")
    print(f"Low quality filtered: {skipped_low_quality}")
    print(f"Number of haplotype patterns: {len(pattern_counts)}")

    # ADR analysis results summary
    print(f"\nADR Analysis Results (threshold={adr_threshold}):")
    print(f"  Genotype: {adr_results['genotype']}")
    print(f"  ADR value: {adr_results['adr_value']:.4f}")
    print(f"  Note: {adr_results['note']}")

    if adr_results['genotype'] == 'heterozygous':
        print(f"  Haplotype1: {adr_results['haplotype1']} (count: {adr_results['count1']})")
        print(f"  Haplotype2: {adr_results['haplotype2']} (count: {adr_results['count2']})")
    elif adr_results['genotype'] == 'homozygous':
        print(f"  Haplotype: {adr_results['haplotype1']}")
        print(f"  Two identical haplotypes, each count: {adr_results['count1']:.1f}")

    print("\nTop 5 most common haplotype patterns:")
    if pattern_counts:
        for i, (pattern, count) in enumerate(pattern_counts.most_common(5), 1):
            percentage = count / covering_reads * 100 if covering_reads > 0 else 0
            print(f"  {i}. {pattern}: {count} reads ({percentage:.2f}%)")

    print("\nAdditional statistics:")
    print(f"  Reads aligned to other references: {other_ref_reads}")
    print(f"  Unmapped reads: {unmapped_reads}")


def main():
    """
    Main function to run SAM file analysis.

    Configure file paths and parameters here.
    """
    # File path configuration
    sam_file = r"D:\MH-001.sam"
    excel_file = r"D:\input_snp_positions.xlsx"
    output_file = r"D:\output_01_haplotype_ranking.xlsx"

    # ADR output file (optional)
    adr_output_file = r"D:\output_02_genotype_call.xlsx"

    # ADR threshold (configurable parameter)
    adr_threshold = 0.1

    print("=" * 60)
    print("SAM File SNP Analysis Tool with ADR Calculation")
    print("=" * 60)
    print(f"Input SAM file: {sam_file}")
    print(f"Input Excel file: {excel_file}")
    print(f"Main output file: {output_file}")
    print(f"ADR output file: {adr_output_file}")
    print(f"ADR threshold: {adr_threshold}")
    print("=" * 60)

    # Execute analysis
    try:
        analyze_single_sam_file(
            sam_file=sam_file,
            excel_file=excel_file,
            output_file=output_file,
            adr_output_file=adr_output_file,
            min_mapq=20,
            adr_threshold=adr_threshold
        )
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
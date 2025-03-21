#!/usr/bin/env python3
import gzip
import re
import argparse
import sys

def parse_mpileup_line(line):
    """Parse a single line from a mpileup file and calculate nucleotide statistics."""
    try:
        fields = line.strip().split('\t')
        
        # Need at least 6 fields (chrom, pos, ref, depth, bases, qualities)
        if len(fields) < 6:
            return None
        
        pos = int(fields[1])
        ref_base = fields[2].upper()
        coverage = int(fields[3])
        read_bases = fields[4]
        
        # If no coverage, return zeros
        if coverage == 0:
            return create_zero_stats(pos, ref_base)
        
        # Clean up the read_bases string
        read_bases = re.sub(r'\^.', '', read_bases)  # Remove read start markers and mapping quality
        read_bases = read_bases.replace('$', '')     # Remove read end markers
        
        # Initialize counters
        matches = 0
        mismatches = 0
        insertions = 0
        deletions = 0
        base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
        
        # Parse the read bases string character by character
        i = 0
        while i < len(read_bases):
            if read_bases[i] in '.,':
                # Match to reference
                matches += 1
                base_counts[ref_base] += 1
                i += 1
            elif read_bases[i].upper() in 'ACGT':
                # Mismatch to reference
                mismatches += 1
                base = read_bases[i].upper()
                base_counts[base] += 1
                i += 1
            elif read_bases[i] == '+':
                # Insertion
                i += 1
                # Extract insertion length
                ins_len_str = ''
                while i < len(read_bases) and read_bases[i].isdigit():
                    ins_len_str += read_bases[i]
                    i += 1
                ins_len = int(ins_len_str) if ins_len_str else 0
                i += ins_len  # Skip the inserted bases
                insertions += 1
            elif read_bases[i] == '-':
                # Deletion
                i += 1
                # Extract deletion length
                del_len_str = ''
                while i < len(read_bases) and read_bases[i].isdigit():
                    del_len_str += read_bases[i]
                    i += 1
                del_len = int(del_len_str) if del_len_str else 0
                i += del_len  # Skip the deleted bases
                deletions += 1
            elif read_bases[i] == '*':
                # Placeholder for a deleted base
                deletions += 1
                i += 1
            else:
                # Skip any other characters
                i += 1
        
        # Calculate totals for events
        total_events = matches + mismatches + insertions + deletions
        
        if total_events == 0:
            return create_zero_stats(pos, ref_base)
        
        # Calculate insertion/deletion percentages relative to total events
        if total_events > 0:
            insertion_percent = round((insertions / total_events) * 100, 4)
            deletion_percent = round((deletions / total_events) * 100, 4)
        else:
            insertion_percent = deletion_percent = 0
            
        # Calculate base percentages
        total_bases = sum(base_counts.values())
        if total_bases > 0:
            # Calculate raw nucleotide percentages
            a_percent_raw = (base_counts['A'] / total_bases) * 100
            g_percent_raw = (base_counts['G'] / total_bases) * 100
            c_percent_raw = (base_counts['C'] / total_bases) * 100
            t_percent_raw = (base_counts['T'] / total_bases) * 100
            
            # Scale the nucleotide percentages so all six percentages sum to 100%
            base_scale_factor = (100 - insertion_percent - deletion_percent) / 100
            
            if base_scale_factor > 0:
                a_percent = round(a_percent_raw * base_scale_factor, 4)
                g_percent = round(g_percent_raw * base_scale_factor, 4)
                c_percent = round(c_percent_raw * base_scale_factor, 4)
                t_percent = round(t_percent_raw * base_scale_factor, 4)
                
                # Fix any rounding errors
                sum_percent = insertion_percent + deletion_percent + a_percent + g_percent + c_percent + t_percent
                if abs(sum_percent - 100) < 0.01:
                    # Find the largest base percentage and adjust it
                    base_percentages = [a_percent, g_percent, c_percent, t_percent]
                    if max(base_percentages) > 0:
                        if a_percent == max(base_percentages):
                            a_percent += (100 - sum_percent)
                        elif g_percent == max(base_percentages):
                            g_percent += (100 - sum_percent)
                        elif c_percent == max(base_percentages):
                            c_percent += (100 - sum_percent)
                        else:
                            t_percent += (100 - sum_percent)
            else:
                # If scale factor is 0, then indels are 100%
                a_percent = g_percent = c_percent = t_percent = 0
            
            # Match percent is the percentage of the reference base
            # Mismatch percent is the second highest nucleotide percentage
            base_dict = {'A': a_percent, 'G': g_percent, 'C': c_percent, 'T': t_percent}
            match_percent = base_dict[ref_base]
            
            # Find the second highest nucleotide percentage (the mismatch)
            other_bases = {b: base_dict[b] for b in base_dict if b != ref_base}
            mismatch_percent = max(other_bases.values()) if other_bases else 0
            
        else:
            a_percent = g_percent = c_percent = t_percent = 0
            match_percent = mismatch_percent = 0
            
            # If no bases but we have indels, make sure they sum to 100%
            if insertion_percent + deletion_percent > 0:
                indel_sum = insertion_percent + deletion_percent
                if abs(indel_sum - 100) < 0.01:
                    if insertion_percent >= deletion_percent:
                        insertion_percent = round(100 - deletion_percent, 4)
                    else:
                        deletion_percent = round(100 - insertion_percent, 4)
        
        return {
            'pos': pos,
            'ref_base': ref_base,
            'coverage': coverage,
            'match_percent': match_percent,
            'mismatch_percent': mismatch_percent,
            'insertion_percent': insertion_percent,
            'deletion_percent': deletion_percent,
            'A_percent': a_percent,
            'G_percent': g_percent,
            'C_percent': c_percent,
            'T_percent': t_percent,
            'raw_bases': matches + mismatches,
            'matches': matches,
            'mismatches': mismatches,
            'insertions': insertions,
            'deletions': deletions,
            'base_counts': base_counts
        }
    except Exception as e:
        # Log the error and return None
        print(f"Error processing line: {line.strip()}", file=sys.stderr)
        print(f"Exception: {str(e)}", file=sys.stderr)
        return None

def create_zero_stats(pos, ref_base):
    """Create a stats dictionary with zero values."""
    return {
        'pos': pos,
        'ref_base': ref_base,
        'coverage': 0,
        'match_percent': 0,
        'mismatch_percent': 0,
        'insertion_percent': 0,
        'deletion_percent': 0,
        'A_percent': 0,
        'G_percent': 0,
        'C_percent': 0,
        'T_percent': 0,
        'raw_bases': 0,
        'matches': 0,
        'mismatches': 0,
        'insertions': 0,
        'deletions': 0,
        'base_counts': {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    }

def process_mpileup(mpileup_file, output_file):
    """Process a gzipped mpileup file and output nucleotide statistics."""
    try:
        # Open the output file
        with open(output_file, 'w') as out_f:
            # Write header
            out_f.write("Ref_Position_1based\tRef_Base\tMatch_Percent\tMismatch_Percent\tInsertion_Percent\tDeletion_Percent\tA_Percent\tG_Percent\tC_Percent\tT_Percent\n")
            
            # Open the gzipped file
            with gzip.open(mpileup_file, 'rt') as f:
                for line_num, line in enumerate(f, 1):
                    try:
                        stats = parse_mpileup_line(line)
                        if stats:
                            # Format all percentage values to exactly 4 decimal places
                            out_f.write(f"{stats['pos']}\t{stats['ref_base']}\t"
                                  f"{stats['match_percent']:.4f}\t"
                                  f"{stats['mismatch_percent']:.4f}\t"
                                  f"{stats['insertion_percent']:.4f}\t"
                                  f"{stats['deletion_percent']:.4f}\t"
                                  f"{stats['A_percent']:.4f}\t"
                                  f"{stats['G_percent']:.4f}\t"
                                  f"{stats['C_percent']:.4f}\t"
                                  f"{stats['T_percent']:.4f}\n")
                    except Exception as e:
                        print(f"Error on line {line_num}: {str(e)}", file=sys.stderr)
                        continue
    except gzip.BadGzipFile:
        print(f"Error: {mpileup_file} is not a valid gzipped file", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File {mpileup_file} not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        sys.exit(1)

def is_polymorphic(stats, threshold=10):
    """Determine if a position is polymorphic (has variation above threshold)."""
    # Position is polymorphic if:
    # 1. Has significant mismatches (above threshold percentage)
    # 2. Has significant insertions or deletions
    return (stats['mismatch_percent'] >= threshold or 
            stats['insertion_percent'] >= threshold or 
            stats['deletion_percent'] >= threshold)

def get_sequence_context(all_stats, pos, window=10):
    """Get the sequence context around a position."""
    ref_seq = ""
    
    # Get reference bases for positions within window
    start_pos = max(1, pos - window)
    end_pos = pos + window
    
    for i in range(start_pos, end_pos + 1):
        if i in all_stats:
            ref_seq += all_stats[i]['ref_base']
        else:
            ref_seq += "N"  # Use N for positions we don't have
    
    # Split the sequence to highlight the reference base
    mid_point = min(pos - start_pos, len(ref_seq) - 1)
    before = ref_seq[:mid_point]
    center = ref_seq[mid_point]
    after = ref_seq[mid_point + 1:]
    
    return before, center, after

def generate_polymorphism_report(mpileup_file, output_file, threshold=10):
    """Generate a report of polymorphic positions."""
    try:
        all_stats = {}
        
        # First pass: read all positions to gather complete stats
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                stats = parse_mpileup_line(line)
                if stats:
                    all_stats[stats['pos']] = stats
        
        # Second pass: identify polymorphic positions and write report
        with open(output_file, 'w') as out_f:
            for pos, stats in sorted(all_stats.items()):
                if is_polymorphic(stats, threshold):
                    # Get sequence context
                    before, center, after = get_sequence_context(all_stats, pos)
                    
                    # Write position and reference information
                    out_f.write(f"(1-based) Position:{pos}, Reference Base={stats['ref_base']}\n")
                    out_f.write(f"Aligned Read Count:{stats['coverage']}\n")
                    
                    # Write table header
                    out_f.write("Mat\tMis\tIns\tDel\tA\tG\tC\tT\n")
                    
                    # Write table values (rounded to whole numbers for readability)
                    out_f.write(f"{round(stats['match_percent'])}\t"
                               f"{round(stats['mismatch_percent'])}\t"
                               f"{round(stats['insertion_percent'])}\t"
                               f"{round(stats['deletion_percent'])}\t"
                               f"{round(stats['A_percent'])}\t"
                               f"{round(stats['G_percent'])}\t"
                               f"{round(stats['C_percent'])}\t"
                               f"{round(stats['T_percent'])}\t\n")
                    
                    # Write sequence context
                    out_f.write(f"{before} {center} {after}\n\n")
                    
    except Exception as e:
        print(f"Error generating polymorphism report: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Analyze mpileup file for nucleotide statistics.')
    parser.add_argument('-i', '--input', required=True, help='Gzipped mpileup file to analyze')
    parser.add_argument('-o', '--output', required=True, help='Output file for nucleotide frequency stats')
    parser.add_argument('-s', '--snps', required=True, help='Output file for SNP summary')
    parser.add_argument('-t', '--threshold', type=int, default=10, 
                        help='Threshold percentage for considering a position polymorphic (default: 10)')
    args = parser.parse_args()
    
    # Process mpileup file to generate nucleotide frequencies
    process_mpileup(args.input, args.output)
    
    # Generate polymorphism report
    generate_polymorphism_report(args.input, args.snps, args.threshold)
    
    print(f"Nucleotide frequency stats written to: {args.output}")
    print(f"Polymorphism report written to: {args.snps}")

if __name__ == "__main__":
    main()
exit(0)

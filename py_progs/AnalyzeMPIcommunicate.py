#!/usr/bin/env python3
"""
MPI Pack/Unpack Sequence Analyzer - Simplified for specific C code format
"""

import re
import sys
from dataclasses import dataclass
from typing import List, Tuple

@dataclass
class MPICall:
    call_type: str  # 'pack' or 'unpack'
    line_number: int
    variable: str
    count: str
    data_type: str

def analyze_mpi_file(filename: str):
    """Analyze MPI Pack/Unpack sequences in the file"""
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    pack_calls = []
    unpack_calls = []
    
    # Extract all MPI calls
    for i, line in enumerate(lines, 1):
        line = line.strip()
        
        if line.startswith('MPI_Pack (') or line.startswith('MPI_Pack('):
            variable, count, data_type = parse_mpi_call(line, 'pack')
            if variable:
                pack_calls.append(MPICall('pack', i, variable, count, data_type))
        
        elif line.startswith('MPI_Unpack (') or line.startswith('MPI_Unpack('):
            variable, count, data_type = parse_mpi_call(line, 'unpack')
            if variable:
                unpack_calls.append(MPICall('unpack', i, variable, count, data_type))
    
    print(f"Found {len(pack_calls)} MPI_Pack calls and {len(unpack_calls)} MPI_Unpack calls\n")
    
    # Check for mismatches
    problems = []
    
    if len(pack_calls) != len(unpack_calls):
        problems.append(f"MAJOR ERROR: {len(pack_calls)} pack calls vs {len(unpack_calls)} unpack calls")
        
    max_compare = min(len(pack_calls), len(unpack_calls))
    
    for i in range(max_compare):
        pack = pack_calls[i]
        unpack = unpack_calls[i]
        
        # Normalize variable names for comparison
        pack_var = normalize_variable(pack.variable)
        unpack_var = normalize_variable(unpack.variable)
        
        # Check for known issues
        if pack.count != unpack.count:
            problems.append(f"Count mismatch at position {i+1}: "
                          f"Pack line {pack.line_number} sends '{pack.count}' but "
                          f"unpack line {unpack.line_number} receives '{unpack.count}'")
        
        if pack.data_type != unpack.data_type:
            problems.append(f"Data type mismatch at position {i+1}: "
                          f"Pack line {pack.line_number} uses '{pack.data_type}' but "
                          f"unpack line {unpack.line_number} uses '{unpack.data_type}'")
        
        # Check for typos in variable names
        if 'ave_freq_norm' in pack.variable and 'ave_freq_nrom' in unpack.variable:
            problems.append(f"TYPO at position {i+1}: "
                          f"Line {unpack.line_number} has 'ave_freq_nrom' (should be 'ave_freq_norm')")
        
        elif pack_var != unpack_var:
            problems.append(f"Variable mismatch at position {i+1}: "
                          f"Pack line {pack.line_number} uses '{pack.variable}' but "
                          f"unpack line {unpack.line_number} uses '{unpack.variable}'")
    
    # Report results
    print("=" * 80)
    print("MPI PACK/UNPACK ANALYSIS RESULTS")
    print("=" * 80)
    
    if not problems:
        print("SUCCESS: All MPI Pack/Unpack sequences match correctly!")
    else:
        print(f"FOUND {len(problems)} PROBLEMS:\n")
        for i, problem in enumerate(problems, 1):
            print(f"{i}. {problem}")
        
        print(f"\n" + "=" * 80)
        print("SPECIFIC ISSUES TO FIX:")
        
        for problem in problems:
            if "TYPO" in problem:
                print(f"- Fix typo: {problem}")
            elif "Count mismatch" in problem:
                print(f"- Fix count: {problem}")
            elif "Data type mismatch" in problem:
                print(f"- Fix data type: {problem}")
            elif "Variable mismatch" in problem:
                print(f"- Fix variable: {problem}")
    
    print("=" * 80)

def parse_mpi_call(line: str, call_type: str) -> Tuple[str, str, str]:
    """Parse an MPI Pack or Unpack call"""
    try:
        # Remove MPI_Pack( or MPI_Unpack( and trailing );
        content = line.replace(f'MPI_{call_type.title()}', '').strip()
        content = content.strip('(').strip(');').strip()
        
        # Split by commas but be careful with nested parentheses
        parts = []
        current = ""
        paren_count = 0
        
        for char in content:
            if char == '(':
                paren_count += 1
            elif char == ')':
                paren_count -= 1
            elif char == ',' and paren_count == 0:
                parts.append(current.strip())
                current = ""
                continue
            current += char
        
        if current.strip():
            parts.append(current.strip())
        
        if call_type == 'pack':
            # MPI_Pack(buffer, count, datatype, packbuf, packbufsize, position, comm)
            if len(parts) >= 3:
                variable = parts[0]    # buffer
                count = parts[1]       # count  
                data_type = parts[2]   # datatype
                return variable, count, data_type
        else:  # unpack
            # MPI_Unpack(packbuf, packbufsize, position, buffer, count, datatype, comm)
            if len(parts) >= 6:
                variable = parts[3]    # buffer
                count = parts[4]       # count
                data_type = parts[5]   # datatype
                return variable, count, data_type
        
    except Exception as e:
        print(f"Error parsing line: {line}")
        print(f"Error: {e}")
    
    return "", "", ""

def normalize_variable(var: str) -> str:
    """Normalize variable for comparison"""
    if not var:
        return var
    
    # Remove address operators
    var = var.replace('&', '').replace('*', '')
    
    # Handle cell-> vs plasmamain[n_plasma]. patterns
    if 'cell->' in var:
        return var.split('->')[-1]
    elif 'plasmamain[' in var and '].' in var:
        return var.split('].')[-1]
    
    return var

def main():
    if len(sys.argv) != 2:
        print("Usage: python mpi_analyzer.py <c_file>")
        sys.exit(1)
    
    try:
        analyze_mpi_file(sys.argv[1])
    except FileNotFoundError:
        print(f"Error: File '{sys.argv[1]}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error analyzing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

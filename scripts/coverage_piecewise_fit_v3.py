#!/usr/bin/env python3.9
# coverage_piecewise_fit_v2.py  written by Bernardo w/ DeepSeek help.   2oct2025    v. 3may2026
# Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno
# to analyse the coverage data of sats

"""
Piecewise Linear Fit for Genomic Read Coverage Data
Fits a two-segment model (inclined + flat) to genomic read coverage data
"""

# usage: coverage_piecewise_fit_v3.py --data_suffix _LILAP_HiFi_   --sat_datafile ../AAAATAT_target3.sat  --min_consecutive_zeros 10 > ../AAAATAT_target3.sat.penet_data



import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import re


def validate_pattern(value):
    if value in ['empirical_0', 'skip_sat','full_ctg']:
        return value
    if re.match(r'^\d+,\d+$', value):
        return value
    raise argparse.ArgumentTypeError(f"Invalid choice: {value}. Allowed: 'empirical_0', 'skip_sat','full_ctg', or 'int,int'")


import pandas as pd

def load_sat_data(filename, delimiter=None):
    """
    Reads a tabular text file with a header row, skipping comment lines starting with '#'.
    Parameters:
        filename (str): Path to the data file.
        delimiter (str, optional): If provided, use this exact delimiter.
                                   If None, automatically handle any whitespace
                                   (tabs, spaces, or combinations) as delimiter.
    Returns:
        pd.DataFrame: DataFrame with column names from the first non‑comment row.
    """
    print('entering load_sat_data function ', filename, file=log_file)
    
    if delimiter is None:
        # delim_whitespace=True treats any whitespace as delimiter
        # comment='#' skips lines that start with '#'
        sat_data = pd.read_csv(filename, delim_whitespace=True, header=0, comment='#')
    else:
        sat_data = pd.read_csv(filename, sep=delimiter, header=0, comment='#')
    
    if 'sequenceID' not in sat_data.columns:
        print("WARNING! file ", filename, " probably has no header. Columns sequenceID, monomer, sat_start, sat_end, and seq_size are required.") 
    return sat_data

def load_depth_data(filename):
    """Load X,Y data from plain text file"""
    try:
        data = np.loadtxt(filename)
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError("File must contain exactly two columns (X and Y)")
        return data[:, 0], data[:, 1]
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")
        sys.exit(1)

def parse_sat_position(sat_position_str):
    """Parse the sat_position string into start and end integers"""
    try:
        start, end = map(int, sat_position_str.split(','))
        if start >= end:
            raise ValueError("Start position must be less than end position")
        return start, end
    except Exception as e:
        print(f"Error parsing sat_position: {e}")
        print("sat_position must be two integers separated by comma (e.g., '1000,2000')")
        sys.exit(1)


def find_first_zero_coverage(x_data, y_data, ref_point, sat_LR, consecutive_zeros):
    """
    Find the first position where there are at least `consecutive` zeros in a row.
    """
    print('########################### entered find_first_zero_coverage', file=log_file)
    # Find all positions with zero coverage
    zero_mask = y_data == 0
    if not np.any(zero_mask):
        return None, None, None    
    zero_positions = x_data[zero_mask]
    count = 0
    
    if sat_LR == "right":
        for i in range(len(y_data)):
            # print('pos:', i, '  cov:' , y_data[i], '  count=', count)
            if y_data[i] == 0:
                count += 1
                if count >= consecutive_zeros:
                    # Found a run of at least 'consecutive' zeros
                    empirical_first_zero_pos = x_data[i - count + 1]   # position of first zero in the run
                    signed_distance = empirical_first_zero_pos - ref_point
                    print('zero block found at pos ', i,file=log_file)              
                    return empirical_first_zero_pos, signed_distance, len(zero_positions)
            else:
                count = 0  # reset counter when a non-zero value appears
            
    if sat_LR == "left":
        for i in range(len(y_data) - 1, -1, -1):
            # print('pos:', i, '  cov:' , y_data[i], '  count=', count)
            if y_data[i] == 0:
                count += 1
                if count >= consecutive_zeros:
                    # Found a run of at least 'consecutive' zeros
                    empirical_first_zero_pos = x_data[i - count + 1]   # position of first zero in the run
                    signed_distance = empirical_first_zero_pos - ref_point
                    print('zero block found at pos ', i, file=log_file)              
                    return empirical_first_zero_pos, signed_distance, len(zero_positions)
            else:
                count = 0  # reset counter when a non-zero value appears     
    
    return None, None, None


def piecewise_quadratic_flat(params, x, y):
    """Objective function for quadratic + flat piecewise fit"""
    breakpoint, a, b, c, flat_level = params
    
    # Ensure breakpoint is within reasonable bounds
    if breakpoint < x.min() or breakpoint > x.max():
        return np.inf
    
    # Calculate predictions
    mask = x < breakpoint
    y_pred = np.where(mask, a*x**2 + b*x + c, flat_level)
    
    # Calculate residual sum of squares
    rss = np.sum((y - y_pred) ** 2)
    return rss

def piecewise_linear_flat(params, x, y):
    """Objective function for linear + flat piecewise fit"""
    breakpoint, slope, intercept, flat_level = params
    
    # Ensure breakpoint is within reasonable bounds
    if breakpoint < x.min() or breakpoint > x.max():
        return np.inf
    
    # Calculate predictions
    mask = x < breakpoint
    y_pred = np.where(mask, slope*x + intercept, flat_level)
    
    # Calculate residual sum of squares
    rss = np.sum((y - y_pred) ** 2)
    return rss

def fit_piecewise_scipy(x_data, y_data, breakpoint=None, quadratic=False, breakpoint_number=1):
    """
    Fit piecewise function using scipy.optimize.minimize
    Returns segments with automatic slope detection
    """
    try:
        from scipy.optimize import minimize
    except ImportError:
        print("Error: scipy library not installed.")
        print("Install with: pip install scipy")
        sys.exit(1)
    
    # For scipy optimization, we currently only support 1 breakpoint
    if breakpoint_number != 1:
        print("Warning: scipy optimization currently only supports 1 breakpoint. Using 1 breakpoint.")
        breakpoint_number = 1
    
    if breakpoint is not None:
        print(f"Using user-specified breakpoint: {breakpoint}", file=log_file)
        # For fixed breakpoint, we only optimize the other parameters
        if quadratic:
            # Initial guesses for quadratic parameters
            a_guess, b_guess, c_guess = np.polyfit(x_data, y_data, 2)
            flat_guess = np.mean(y_data)
            initial_params = [a_guess, b_guess, c_guess, flat_guess]
            
            # Define objective function for fixed breakpoint
            def objective_fixed(params):
                a, b, c, flat_level = params
                mask = x_data < breakpoint
                y_pred = np.where(mask, a*x_data**2 + b*x_data + c, flat_level)
                return np.sum((y_data - y_pred) ** 2)
            
            # Optimize
            result = minimize(objective_fixed, initial_params, method='L-BFGS-B')
            a_opt, b_opt, c_opt, flat_opt = result.x
            slopes = [0, 0]  # Not used for quadratic
            intercepts = [c_opt, flat_opt]
            
        else:
            # Initial guesses for linear parameters
            slope_guess, intercept_guess = np.polyfit(x_data, y_data, 1)
            flat_guess = np.mean(y_data)
            initial_params = [slope_guess, intercept_guess, flat_guess]
            
            # Define objective function for fixed breakpoint
            def objective_fixed(params):
                slope, intercept, flat_level = params
                mask = x_data < breakpoint
                y_pred = np.where(mask, slope*x_data + intercept, flat_level)
                return np.sum((y_data - y_pred) ** 2)
            
            # Optimize
            result = minimize(objective_fixed, initial_params, method='L-BFGS-B')
            slope_opt, intercept_opt, flat_opt = result.x
            slopes = [slope_opt, 0]
            intercepts = [intercept_opt, flat_opt]
            
        final_breakpoint = breakpoint
        
    else:
        # Optimize breakpoint and all parameters
        if quadratic:
            # Initial guesses
            breakpoint_guess = np.median(x_data)
            a_guess, b_guess, c_guess = np.polyfit(x_data, y_data, 2)
            flat_guess = np.mean(y_data)
            initial_params = [breakpoint_guess, a_guess, b_guess, c_guess, flat_guess]
            
            # Bounds for parameters
            bounds = [
                (x_data.min(), x_data.max()),  # breakpoint
                (-np.inf, np.inf),  # a
                (-np.inf, np.inf),  # b
                (-np.inf, np.inf),  # c
                (y_data.min(), y_data.max())  # flat_level
            ]
            
            # Optimize
            result = minimize(piecewise_quadratic_flat, initial_params, args=(x_data, y_data),
                            method='L-BFGS-B', bounds=bounds)
            final_breakpoint, a_opt, b_opt, c_opt, flat_opt = result.x
            slopes = [0, 0]  # Not used for quadratic
            intercepts = [c_opt, flat_opt]
            
        else:
            # Initial guesses
            breakpoint_guess = np.median(x_data)
            slope_guess, intercept_guess = np.polyfit(x_data, y_data, 1)
            flat_guess = np.mean(y_data)
            initial_params = [breakpoint_guess, slope_guess, intercept_guess, flat_guess]
            
            # Bounds for parameters
            bounds = [
                (x_data.min(), x_data.max()),  # breakpoint
                (-np.inf, np.inf),  # slope
                (-np.inf, np.inf),  # intercept
                (y_data.min(), y_data.max())  # flat_level
            ]
            
            # Optimize
            result = minimize(piecewise_linear_flat, initial_params, args=(x_data, y_data),
                            method='L-BFGS-B', bounds=bounds)
            final_breakpoint, slope_opt, intercept_opt, flat_opt = result.x
            slopes = [slope_opt, 0]
            intercepts = [intercept_opt, flat_opt]
    
    # Calculate predictions and R-squared
    if quadratic:
        mask = x_data < final_breakpoint
        y_pred = np.where(mask, a_opt*x_data**2 + b_opt*x_data + c_opt, flat_opt)
    else:
        mask = x_data < final_breakpoint
        y_pred = np.where(mask, slope_opt*x_data + intercept_opt, flat_opt)
    
    ss_res = np.sum((y_data - y_pred) ** 2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
    
    # For scipy optimization, we always assume inclined then flat
    inclined_segment = 0
    flat_segment = 1
    
    result_dict = {
        'breakpoints': [final_breakpoint],
        'slopes': slopes,
        'intercepts': intercepts,
        'segment_types': ['inclined', 'flat'],
        'r_squared': r_squared,
        'predictions': y_pred,
        'quadratic': quadratic
    }
    
    if quadratic:
        result_dict.update({
            'quadratic_coeffs': [a_opt, b_opt, c_opt],
            'quadratic_poly': np.poly1d([a_opt, b_opt, c_opt]),
            'flat_level': flat_opt
        })
    else:
        result_dict.update({
            'flat_level': flat_opt
        })
    
    return result_dict

def fit_piecewise_pwlf(x_data, y_data, breakpoint=None, quadratic=False, breakpoint_number=1):
    """
    Fit piecewise function using pwlf library with multiple breakpoints
    Returns segments with automatic slope detection
    """
    try:
        import pwlf
    except ImportError:
        print("Error: pwlf library not installed.")
        print("Install with: pip install pwlf")
        sys.exit(1)
    
    # Initialize
    my_pwlf = pwlf.PiecewiseLinFit(x_data, y_data)
    
    # Calculate number of segments (breakpoints + 1)
    n_segments = breakpoint_number + 1
    
    if breakpoint is not None:
        # Use user-specified breakpoint(s)
        if isinstance(breakpoint, (int, float)):
            # Single breakpoint provided
            print(f"Using user-specified breakpoint: {breakpoint}", file=log_file)
            custom_breaks = [x_data.min(), breakpoint, x_data.max()]
        else:
            # Multiple breakpoints provided
            print(f"Using user-specified breakpoints: {breakpoint}", file=log_file)
            custom_breaks = [x_data.min()] + sorted(breakpoint) + [x_data.max()]
        
        breaks = my_pwlf.fit_with_breaks(custom_breaks)
    else:
        # Let pwlf find the optimal breakpoints
        print(f"Finding {breakpoint_number} breakpoints ({n_segments} segments)...", file=log_file)
        breaks = my_pwlf.fit(n_segments)
    
    # Get slopes and intercepts
    slopes = my_pwlf.calc_slopes()
    
    # Try different ways to get intercepts based on pwlf version
    try:
        # Newer versions: intercepts is a property
        intercepts = my_pwlf.intercepts
    except TypeError:
        # Older versions: intercepts might be a method
        intercepts = my_pwlf.intercepts()
    
    # If quadratic is requested, we can only handle 1 breakpoint
    if quadratic:
        if breakpoint_number > 1:
            print("Warning: Quadratic fitting currently only supports 1 breakpoint. Using 1 breakpoint.")
            breakpoint_number = 1
        
        print("Fitting quadratic curve to inclined segment...", file=log_file)
        # For quadratic, we assume the first segment is inclined and the rest are flat
        # Get data for the first segment
        inclined_mask = x_data < breaks[1]
        x_inclined = x_data[inclined_mask]
        y_inclined = y_data[inclined_mask]
        
        # Fit quadratic curve to the first segment
        quadratic_coeffs = np.polyfit(x_inclined, y_inclined, 2)
        quadratic_poly = np.poly1d(quadratic_coeffs)
        
        # Calculate predictions for inclined segment
        y_pred_inclined = quadratic_poly(x_inclined)
        
        # For other segments, use the original linear fit
        y_pred = my_pwlf.predict(x_data)
        y_pred[inclined_mask] = y_pred_inclined
        
        # Update intercepts for display
        flat_level = np.mean(y_data[~inclined_mask])
        intercepts = [intercepts[0], flat_level]
        
        # Calculate R-squared for the combined fit
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # For quadratic with 1 breakpoint, we have inclined then flat
        segment_types = ['quadratic', 'flat']
        
        return {
            'breakpoints': breaks[1:-1],
            'slopes': slopes,
            'intercepts': intercepts,
            'segment_types': segment_types,
            'r_squared': r_squared,
            'predictions': y_pred,
            'model': my_pwlf,
            'quadratic': True,
            'quadratic_coeffs': quadratic_coeffs,
            'quadratic_poly': quadratic_poly,
            'flat_level': flat_level
        }
    else:
        # Standard linear fit with automatic slope detection
        y_pred = my_pwlf.predict(x_data)
        ss_res = np.sum((y_data - y_pred) ** 2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
        
        # Automatically detect segment types based on slopes
        slope_threshold = 0.001
        segment_types = []
        
        for i, slope in enumerate(slopes):
            if abs(slope) <= slope_threshold:
                segment_types.append('flat')
            else:
                segment_types.append('inclined')
        
        return {
            'breakpoints': breaks[1:-1],
            'slopes': slopes,
            'intercepts': intercepts,
            'segment_types': segment_types,
            'r_squared': r_squared,
            'predictions': y_pred,
            'model': my_pwlf,
            'quadratic': False
        }

def calculate_x_intercept_linear(slope, intercept):
    """Calculate X-axis intercept (where y=0) for linear segment"""
    if abs(slope) < 1e-10:
        return None
    return -intercept / slope

def calculate_x_intercept_quadratic(quadratic_coeffs):
    """Calculate X-axis intercepts (where y=0) for quadratic curve"""
    a, b, c = quadratic_coeffs
    
    # Solve quadratic equation: ax² + bx + c = 0
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        # No real roots
        return None, None
    
    sqrt_discriminant = np.sqrt(discriminant)
    x1 = (-b - sqrt_discriminant) / (2*a)
    x2 = (-b + sqrt_discriminant) / (2*a)
    
    return x1, x2

def create_plot(x_data, y_data, result, output_file, format_type, sat_start, sat_end, 
                ref_point=None, x_intercept=None, x_intercept_quadratic=None, 
                breakpoint_ref_distance=None, empirical_first_zero_pos=None, graph_ref_lines=0):
    """Create and save visualization of the fit"""
    plt.figure(figsize=(12, 8))
    
    # plt.ylim(bottom=0) # show y >= 0  ficou mto ruim
    
    # Plot original data
    plt.scatter(x_data, y_data, alpha=0.6, label='Data', s=10, color='blue')
    
    # Plot fitted model
    x_fit = np.linspace(x_data.min(), x_data.max(), 1000)
    
    
    
    
    if result.get('quadratic', False):
        # For quadratic fit, we need to create the prediction manually
        bp = result['breakpoints'][0]  # Only one breakpoint for quadratic
        quadratic_poly = result['quadratic_poly']
        flat_level = result.get('flat_level', result['intercepts'][1])
        
        # Inclined segment first, then flat
        mask_inclined = x_fit < bp
        y_fit = np.where(mask_inclined, quadratic_poly(x_fit), flat_level)
    else:
        # Standard linear fit
        y_fit = result['model'].predict(x_fit)
    
    
    
    # Create mask for non‑negative y values
    mask = y_fit >= 0
    # Plot only the points where y_fit >= 0
    plt.plot(x_fit[mask], y_fit[mask], 'r-', linewidth=2, label='Fitted model')
    
    
    # plt.plot(x_fit, y_fit, 'r-', linewidth=2, label='Fitted model')
    
    # Mark all breakpoints
    breakpoints = result['breakpoints']
    for i, bp in enumerate(breakpoints):
        plt.axvline(x=bp, color='green', linestyle='--', 
                    label=f'Breakpoint {i+1}: {bp:.0f}  [{relative_breakpoint_pos}]')
    
    # Mark first zero coverage position if available
    if empirical_first_zero_pos is not None:
        plt.axvline(x=empirical_first_zero_pos, color='brown', linestyle='dotted', linewidth=2, label=f'First zero: {empirical_first_zero_pos:.0f}  [{relative_empirical_first_zero_pos}]')
    
    # Add reference point and related lines based on graph_ref_lines setting
    if ref_point is not None and graph_ref_lines >= 1:
        # Mark reference point
        plt.axvline(x=ref_point, color='orange', linestyle='-', 
                   linewidth=2, label=f'Reference: {ref_point}')
        
        # Add X-intercept lines if graph_ref_lines >= 2
        if graph_ref_lines >= 2:
            # For linear fits
            if x_intercept is not None and not result.get('quadratic', False):
                # Mark X-intercept
                plt.axvline(x=x_intercept, color='purple', linestyle='-', 
                           linewidth=2, label=f'X-intercept: {x_intercept:.2f}')
                
                # Draw distance line between reference and X-intercept
                y_pos = plt.ylim()[1] * 0.7
                plt.plot([ref_point, x_intercept], [y_pos, y_pos], 
                        color='black', linewidth=2, linestyle=':')
                
                # Add distance annotation
                distance = x_intercept - ref_point  # Signed distance
                mid_point = (ref_point + x_intercept) / 2
                direction = "right" if distance > 0 else "left"
                plt.annotate(f'Distance: {distance:.2f} ({direction})', 
                            xy=(mid_point, y_pos),
                            xytext=(0, 10), textcoords='offset points',
                            ha='center', va='bottom', fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            # For quadratic fits - show only the closest intercept
            elif x_intercept_quadratic is not None and result.get('quadratic', False):
                # Mark only the closest X-intercept for quadratic
                plt.axvline(x=x_intercept_quadratic, color='purple', linestyle='-', 
                           linewidth=2, label=f'X-intercept: {x_intercept_quadratic:.2f}')
                
                # Draw distance line between reference and X-intercept
                y_pos = plt.ylim()[1] * 0.7
                plt.plot([ref_point, x_intercept_quadratic], [y_pos, y_pos], 
                        color='black', linewidth=2, linestyle=':')
                
                # Add distance annotation
                distance = x_intercept_quadratic - ref_point  # Signed distance
                mid_point = (ref_point + x_intercept_quadratic) / 2
                direction = "right" if distance > 0 else "left"
                plt.annotate(f'Distance: {distance:.2f} ({direction})', 
                            xy=(mid_point, y_pos),
                            xytext=(0, 10), textcoords='offset points',
                            ha='center', va='bottom', fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Add SAT position rectangle if provided
    if sat_start is not None:
        y_min, y_max = plt.ylim()
        
        # Calculate rectangle position (at bottom, with height = 2% of y-range)
        rect_height = (y_max - y_min) * 0.02
        rect_bottom = y_min - rect_height * 1.5  # Position below the x-axis
        
        plt.axvspan(sat_start, sat_end, ymin=0, ymax=0.05, 
                   facecolor='red', alpha=0.7, label='SAT Position' + ' (' + str(sat_start) +  '-' + str(sat_end) + ')')
        
        # Adjust y-axis to make room for the rectangle
        plt.ylim(rect_bottom, y_max)
    
    plt.xlabel('Genomic Position')
    plt.ylabel('Read Coverage')
    
    title = f'Piecewise Fit to Genomic Read Coverage ({len(breakpoints)} breakpoints)'
    if result.get('quadratic', False):
        title += ' (Quadratic Inclined)'
    plt.title(title)
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add statistics box
    stats_text = f'Breakpoints: {", ".join([f"{bp:.2f}" for bp in breakpoints])}\n'
    
    # Define segment_types for use in statistics
    segment_types = result.get('segment_types', [])
    
    if result.get('quadratic', False):
        a, b, c = result['quadratic_coeffs']
        flat_level = result.get('flat_level', result['intercepts'][1])
        stats_text += f'Quadratic: {a:.4f}x² + {b:.4f}x + {c:.4f}\n'
        stats_text += f'Flat level: {flat_level:.4f}\n'
    else:
        # Show slopes for each segment
        for i, (slope, seg_type) in enumerate(zip(result['slopes'], segment_types)):
            stats_text += f'{seg_type.capitalize()} slope {i+1}: {slope:.4f}\n'
    
    stats_text += f'R²: {result["r_squared"]:.4f}'
    
    if sat_start is not None:
        stats_text += f'\nSAT: {sat_start}-{sat_end}'
    
    if ref_point is not None:
        if not result.get('quadratic', False) and x_intercept is not None:
            distance = x_intercept - ref_point  # Signed distance
            direction = "right" if distance > 0 else "left"
            stats_text += f'\nRef→Xintercept: {distance:.2f} ({direction})'
        elif result.get('quadratic', False) and x_intercept_quadratic is not None:
            distance = x_intercept_quadratic - ref_point  # Signed distance
            direction = "right" if distance > 0 else "left"
            stats_text += f'\nRef→Xintercept: {distance:.2f} ({direction})'
        
        # Add breakpoint-reference distance (signed) for closest breakpoint
        if breakpoint_ref_distance is not None:
            direction_bp = "right" if breakpoint_ref_distance > 0 else "left"
            stats_text += f'\nRef→Breakpoint: {breakpoint_ref_distance:.2f} ({direction_bp})'
        
        # Add first zero coverage distance (signed)
        if empirical_first_zero_pos is not None:
            zero_distance = empirical_first_zero_pos - ref_point  # Signed distance
            direction_zero = "right" if zero_distance > 0 else "left"
            stats_text += f'\nRef→FirstZero: {zero_distance:.2f} ({direction_zero})'
    
    '''
    plt.annotate(stats_text, xy=(0.02, 0.98), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                verticalalignment='top', fontfamily='monospace')
    '''
    plt.tight_layout()
    
    # Save plot
    if format_type.lower() == 'svg':
        plt.savefig(output_file, format='svg', dpi=300)
    elif format_type.lower() == 'png':
        plt.savefig(output_file, format='png', dpi=300)
    elif format_type.lower() == 'pdf':
        plt.savefig(output_file, format='pdf')
    
    print(f"Plot saved as: {output_file}", file=log_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Fit piecewise function with multiple breakpoints to genomic read coverage data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s data.txt
  %(prog)s data.txt --start 1000 --end 5000
  %(prog)s data.txt --graph png
  %(prog)s data.txt --start 1000 --end 5000 --graph svg
  %(prog)s data.txt --sat_position 1500,2500 --graph png
  %(prog)s data.txt --ref 2000 --graph png
  %(prog)s data.txt --breakpoint 3000 --graph png
  %(prog)s data.txt --quadratic --graph png
  %(prog)s data.txt --graph png --graph_suffix "_test1"
  %(prog)s data.txt --optimization scipy --quadratic --graph png
  %(prog)s data.txt --breakpoint_number 2 --graph png
  %(prog)s data.txt --breakpoint_number 3 --graph_ref_lines 1 --graph png
        """
    )
    
    # parser.add_argument('datafile', help='Input data file (two columns: X Y)')
    parser.add_argument('--depth_datafile', type=str, help='Input data file. Two columns: X (position) Y (coverage)')
    parser.add_argument('--start', type=float, help='Start position for piecewise reg. analysis')
    parser.add_argument('--end', type=float, help='End position for piecewise reg. analysis')
    parser.add_argument('--graph', choices=['none', 'svg', 'png', 'pdf'], 
                       default='png', help='Output graph format (default: none)')
    parser.add_argument('--sat_position', type=str, help='SAT position as two integers separated by comma (e.g., "1500,2500")')
    parser.add_argument('--ref', type=float, help='Reference position to calculate distance to X-intercept of inclined line')
    parser.add_argument('--breakpoint', type=float, help='Use specified breakpoint instead of fitting it')
    parser.add_argument('--quadratic', action='store_true',
                       help='Fit quadratic curve to inclined segment instead of linear')
    parser.add_argument('--graph_suffix', type=str, default='v3',
                       help='Suffix to add to graph filename to avoid overwriting')
    parser.add_argument('--optimization', choices=['pwlf', 'scipy'], default='pwlf',
                       help='Optimization method to use (default: pwlf)')
    parser.add_argument('--breakpoint_number', type=int, default=1, choices=[1, 2, 3],
                       help='Number of breakpoints to fit (default: 1)')
    parser.add_argument('--graph_ref_lines', type=int, default=0, choices=[0, 1, 2],
                       help='Level of reference lines in graph: 0=none, 1=ref only, 2=all (default: 0)')
    parser.add_argument('--min_consecutive_zeros', type=int, default=1, 
                       help='Will detect zero coverage only when x consecutive zeros(to avoid triggering a zero coverage point due to a limited drop in coverage. default  is 1 for backward compatibility)')        
    parser.add_argument('--regn_limits', type=validate_pattern, default='empirical_0',  help='start and end limits for piecewise regression. options: skip_sat , empirical_0, full_ctg,  1,30000 [custom]')
    parser.add_argument('--find0_limits', type=str, default='full_ctg',  help='start and end limits for empit=rical search of read zero coverage . options: full_ctg,  1,30000 [custom]')
    parser.add_argument('--sat_datafile', type=str,  help='tabular file containing satelliete locations')
    parser.add_argument('--no-batch', dest='batch_mode', action='store_false', default=True,
                    help='Disable batch mode (batch mode is enabled by default)')
    parser.add_argument('--data_suffix', type=str, default='_ONT_HiFi_',  help='for use in bacth mode. adapts the name of the depth file')
  
    global args
    args = parser.parse_args()
    
       
    
    # aqui entra  o  modo novo, de ler sat_file
    if args.batch_mode:
        # Check if sat file exists
        if not os.path.exists(args.sat_datafile):
            print(f"Error: File '{args.sat_datafile}' not found")
            sys.exit(1)
        else:
            print(f"{'sequenceID':<30}{'monomer':>10}{'sat_LR':>10}{'ref_point':>10}{'asm_penet':>13}{'read_penet':>13}{'x_intercept':>13}{'breakpt':>10}")
            # Open the log file
            global log_file
            log_file = open("coverage_piecewise_v3.log", "w")   # or 'a' for append
 
            sat_data = load_sat_data(args.sat_datafile) # needs the header!!!!!!
            # sequenceID                             monomer   mon_bp   copies    start      end     size   orient   seq_bp   m_blks sat_position
            # ptg000005l:255001-299723                 AATAT        5     3320    17439    44714    27276        R    44723      262
            # print(sat_data)
            for idx, row in sat_data.iterrows():
                # row is a Series, column names are keys
                global sequenceID, monomer, sat_start,sat_end , seq_size
                sequenceID = row.sequenceID
                monomer = row.monomer
                sat_start = row.start
                sat_end = row.end
                seq_size  = row.seq_bp
                #adjusting the name for the depth file
                # ptg000041l:1-50080 become  ptg000041l_0_50080_LILAP_HiFi_AATATAT  # note 0-based issue
                s = sequenceID 
                sID_edited = f"{s.split(':')[0]}_{int(s.split(':')[1].split('-')[0]) - 1}_{s.split(':')[1].split('-')[1]}"
                depth_file = sID_edited + args.data_suffix + monomer + ".depth"   
                # find if sat block is on the left or right of the main sequence.
                global sat_LR
                sat_middle = (sat_start + sat_end)/2
                if  seq_size/2  > (sat_start + sat_end)/2:
                    sat_LR = "left"		
                    ref_point = sat_end			
                    if args.regn_limits == "skip_sat":
                        regress_start = sat_end  
                        regress_end = seq_size
                    if args.regn_limits == "full_ctg":
                        regress_start = 1  
                        regress_end = seq_size  
                    if re.match(r'^\d+,\d+$', args.regn_limits):   #  exact coordinates: 2000,10000
                        parts = args.regn_limits.split(',')
                        regress_start = int(parts[0])
                        regress_end = int(parts[1])
                    if args.regn_limits == "empirical_0":  # this will be dealth with inside the execute function
                        regress_start = 0  
                        regress_end = 0                 
                else:
                    sat_LR = "right"
                    ref_point = sat_start	
                    if args.regn_limits == "skip_sat":
                        regress_start = 1  
                        regress_end = sat_start
                    if args.regn_limits == "full_ctg":
                        regress_start = 1  
                        regress_end = seq_size  
                    if re.match(r'^\d+,\d+$', args.regn_limits):   #  exact coordinates: 2000,10000
                        parts = args.regn_limits.split(',')
                        regress_start = int(parts[0])
                        regress_end = int(parts[1])
                    if args.regn_limits == "empirical_0":  # this will be dealth with inside the execute function
                        regress_start = 0  
                        regress_end = 0                 		

                execute(depth_file,regress_start,regress_end,
                        ref_point, sat_start,sat_end)   
            log_file.close()         
    else:
        # aqui entra o modo antigo, com os patrametros antigos, para usar acoplado ao get_sat_stats_v2.awk 
        # Parse SAT position if provided
        sat_position = None
        if args.sat_position:
            sat_position = parse_sat_position(args.sat_position)
            print(f"SAT position: {sat_position[0]} - {sat_position[1]}", file=log_file)
        else:
            print(f"Error: parameter --sat_position required for non-batch mode.")
            sys.exit(1)
        execute(args.depth_datafile,args.start,args.end,
                args.ref, sat_position[0],sat_position[1])
    
      
def execute(depth_data,regress_start=None,regress_end=None, ref_point=None, sat_start=None,sat_end=None):
    # othre parameters  will be passed by the global args
    # coverage_piecewise_fit_v3temp.py  ptg000005l_255000_299723_ONT_HiFi_AATAT.depth  --start 1 --end 18439 --ref 17439 --sat_position 17439,44714 --graph png --graph_suffix v2_pwlf --optimization pwlf --breakpoint_number 1 --graph_ref_lines 0 --min_consecutive_zeros 10
    
    # Load depth  data
    print(f"Loading data from {depth_data}...", file=log_file)
    x_full, y_full = load_depth_data(depth_data)
    

    # Find first zero coverage position if ref_point is provided
    empirical_first_zero_pos = None
    first_zero_distance = None
    n_zero_points = 0
    
    if ref_point is not None:
        # empirical_first_zero_pos, first_zero_distance, n_zero_points = find_first_zero_coverage(x_full_data, y_full_data, ref_point,sat_LR,args.min_consecutive_zeros)
        empirical_first_zero_pos, first_zero_distance, n_zero_points = find_first_zero_coverage(x_full, y_full, ref_point,sat_LR,args.min_consecutive_zeros)

        if empirical_first_zero_pos is not None:
            print(f"First zero coverage position: {empirical_first_zero_pos:.2f} ({n_zero_points} zero-coverage points found)", file=log_file)
        else:
            print("No zero coverage points found in the dataset", file=log_file)
    
    if args.regn_limits == "empirical_0" and sat_LR == 'left':
        if empirical_first_zero_pos is not None:
            regress_start = empirical_first_zero_pos
        else:
            regress_start = 1    
        regress_end = seq_size      
    if args.regn_limits == "empirical_0" and sat_LR == 'right':
        regress_start = 1
        if empirical_first_zero_pos is not None:
            regress_end = empirical_first_zero_pos       
        else:   
            regress_end = seq_size    

    # Filter data based on start/end parameters
    if regress_start is not None or regress_end is not None:
        mask = np.ones_like(x_full, dtype=bool)
        if regress_start is not None:
            mask = mask & (x_full >= regress_start)
        if regress_end is not None:
            mask = mask & (x_full <= regress_end)
        
        if np.sum(mask) == 0:
            print("Error: No data points in the specified range")
            sys.exit(1)
            
        x_data = x_full[mask]
        y_data = y_full[mask]
        print(f"Using {len(x_data)} points from range {x_data.min():.2f} to {x_data.max():.2f}", file=log_file)
    else:
        x_data = x_full
        y_data = y_full
        print(f"Using all {len(x_data)} data points", file=log_file)
    
    # Check if we have enough data
    min_points_per_segment = 4
    min_total_points = min_points_per_segment * (args.breakpoint_number + 1)
    # print('')
    if len(x_data) < min_total_points:
        print(f"Error: Need at least {min_total_points} data points for {args.breakpoint_number} breakpoints")
        sys.exit(1)
    
    # Sort data by x (required by pwlf)
    sort_idx = np.argsort(x_data)
    x_data = x_data[sort_idx]
    y_data = y_data[sort_idx]
    
    # full data 
    sort_idx = np.argsort(x_full)
    x_full_data = x_full[sort_idx]
    y_full_data = y_full[sort_idx]   
    
        
    
    # Perform the fit
    print("Fitting piecewise model...", file=log_file)
    if args.optimization == 'scipy':
        print("Using scipy direct optimization method...", file=log_file)
        result = fit_piecewise_scipy(x_data, y_data, breakpoint=args.breakpoint, 
                                   quadratic=args.quadratic, breakpoint_number=args.breakpoint_number)
    else:
        print("Using pwlf optimization method...", file=log_file)
        result = fit_piecewise_pwlf(x_data, y_data, breakpoint=args.breakpoint, 
                                  quadratic=args.quadratic, breakpoint_number=args.breakpoint_number)
    
    # Find breakpoint closest to ref_point (if provided)
    closest_breakpoint = None
    breakpoint_ref_distance = None
    
    # Check if breakpoints array is not empty (fix for numpy array ambiguity)
    if ref_point is not None and len(result['breakpoints']) > 0:
        # Find the breakpoint closest to the ref_point
        closest_breakpoint = min(result['breakpoints'], key=lambda bp: abs(bp - ref_point))
        breakpoint_ref_distance = closest_breakpoint - ref_point  # Signed distance
        print(f"Breakpoint closest to ref_point: {closest_breakpoint:.4f}", file=log_file)
    
    # Calculate X-intercept and distance to ref_point if --ref provided
    x_intercept = None
    x_intercept_quadratic = None
    distance = None
    
    if ref_point is not None:
        if args.quadratic:
            # Calculate X-intercept for quadratic curve
            quadratic_coeffs = result['quadratic_coeffs']
            x_int1, x_int2 = calculate_x_intercept_quadratic(quadratic_coeffs)
            # Find the closest intercept to the ref_point point
            if x_int1 is not None or x_int2 is not None:
                intercepts = [xi for xi in [x_int1, x_int2] if xi is not None]
                if intercepts:
                    # Find the closest intercept that's within a reasonable range
                    data_range = (x_data.min(), x_data.max())
                    range_extension = (data_range[1] - data_range[0]) * 0.5  # 50% extension
                    extended_range = (data_range[0] - range_extension, data_range[1] + range_extension)
                    # Filter intercepts to those within a reasonable range
                    valid_intercepts = [xi for xi in intercepts if extended_range[0] <= xi <= extended_range[1]]
                    if valid_intercepts:
                        # Use the closest valid intercept
                        x_intercept_quadratic = min(valid_intercepts, key=lambda xi: abs(xi - ref_point))
                        distance = x_intercept_quadratic - ref_point  # Signed distance
                    else:
                        # If no intercepts in reasonable range, use the closest one
                        x_intercept_quadratic = min(intercepts, key=lambda xi: abs(xi - ref_point))
                        distance = x_intercept_quadratic - ref_point  # Signed distance
                        print("Warning: X-intercept is outside the reasonable range of the data")
                else:
                    print("Warning: No real X-intercepts found for quadratic curve")
            else:
                print("Warning: No real X-intercepts found for quadratic curve")
        else:
            # For multiple segments, find the first inclined segment for X-intercept calculation
            segment_types = result.get('segment_types', [])
            slopes = result['slopes']
            intercepts = result['intercepts']
            
            # Find the first inclined segment
            inclined_idx = None
            for i, seg_type in enumerate(segment_types):
                if seg_type == 'inclined':
                    inclined_idx = i
                    break
            
            if inclined_idx is not None:
                slope = slopes[inclined_idx]
                intercept = intercepts[inclined_idx]
                
                x_intercept = calculate_x_intercept_linear(slope, intercept)
                
                if x_intercept is not None:
                    distance = x_intercept - ref_point  # Signed distance
                else:
                    print("Warning: Could not calculate X-intercept (inclined line is nearly horizontal)")
            else:
                print("No inclined segments found for X-intercept calculation", file=log_file)
    
    # Print results
    print("\n" + "="*50, file=log_file)
    print("PIECEWISE FIT RESULTS", file=log_file)
    print(f"Optimization method: {args.optimization}", file=log_file)
    print(f"Number of breakpoints: {args.breakpoint_number}", file=log_file)
    if args.breakpoint:
        print(f"(Breakpoint fixed at: {args.breakpoint})", file=log_file)
    if args.quadratic:
        print("(Quadratic fit for inclined segment)", file=log_file)
    print("="*50, file=log_file)
    
    breakpoints = result['breakpoints']
    segment_types = result.get('segment_types', [])
    
    print(f"Breakpoints (genomic positions): {', '.join([f'{bp:.4f}' for bp in breakpoints])}", file=log_file)
    
    if ref_point is not None and closest_breakpoint is not None:
        print(f"Breakpoint closest to ref_point: {closest_breakpoint:.4f}", file=log_file)
    
    # Print segment information
    print(f"\nSegment information:", file=log_file)
    for i, (seg_type, slope, intercept) in enumerate(zip(segment_types, result['slopes'], result['intercepts'])):
        if args.quadratic and i == 0:
            a, b, c = result['quadratic_coeffs']
            print(f"Segment {i+1} (Quadratic):", file=log_file)
            print(f"  Equation: y = {a:.6f}x² + {b:.6f}x + {c:.6f}")
        else:
            print(f"Segment {i+1} ({seg_type.capitalize()}):", file=log_file)
            print(f"  Equation: y = {slope:.6f}x + {intercept:.6f}", file=log_file)
            print(f"  Slope: {slope:.6f}", file=log_file)
    
    print(f"\nGoodness of fit:", file=log_file)
    print(f"  R-squared: {result['r_squared']:.6f}", file=log_file)
    print(f"  Explained variance: {result['r_squared']*100:.2f}%", file=log_file)
    
    if sat_start  is not None:
        print(f"\nSAT Position: {sat_start} - {sat_end}", file=log_file)
    
    if ref_point is not None:
        print(f"\nREFERENCE POINT ANALYSIS:", file=log_file)
        print(f"  ref_point position: {ref_point}", file=log_file)
        
        if closest_breakpoint is not None:
            print(f"  Closest breakpoint position: {closest_breakpoint:.4f}", file=log_file)
            print(f"  Breakpoint - ref_point: {breakpoint_ref_distance:.4f} ({'right' if breakpoint_ref_distance > 0 else 'left'} of ref_point)", file=log_file)
        
        if empirical_first_zero_pos is not None:
            print(f"  First zero coverage position: {empirical_first_zero_pos:.4f}", file=log_file)
            print(f"  First zero coverage - ref_point: {first_zero_distance:.4f} ({'right' if first_zero_distance > 0 else 'left'} of ref_point)", file=log_file)
            print(f"  Number of zero coverage points: {n_zero_points}", file=log_file)
        else:
            print(f"  No zero coverage points found in dataset", file=log_file)
        
        if args.quadratic:
            if x_intercept_quadratic is not None:
                print(f"  Closest X-intercept: {x_intercept_quadratic:.4f}", file=log_file)
                print(f"  X-intercept - ref_point: {distance:.4f} ({'right' if distance > 0 else 'left'} of ref_point)", file=log_file)
            else:
                print("  No suitable X-intercept found for quadratic curve", file=log_file)
        else:
            if x_intercept is not None:
                print(f"  X-intercept of first inclined segment: {x_intercept:.4f}", file=log_file)
                print(f"  X-intercept - ref_point: {distance:.4f} ({'right' if distance > 0 else 'left'} of ref_point)", file=log_file)
            else:
                print("  Could not calculate X-intercept (no inclined segments or line is nearly horizontal)", file=log_file)
    
      

 
    # Generate main output data 
    # print(f"{'sequenceID':<30}{'monomer':>10}{'sat_LR':>10}{'ref_point':>10}"
      # f"{'asm_penet':>13}{'read_penet':>12}{'breakpt':>10}")
    
    # in v3 I changed the definition of asm_penet (in v2 was called LILAP_asm)
    # LILAP_asm for left was simply the ref_point ie, the sat_end). For right, (seq_size - ref_point + 1)
    # v3 asm_penet will be simply the sat block size: sat_end - sat_start +1 
    # for a toxic sat in  LILAP the two definitions will be identical or nealry so. |But for the ONT_hifiasm assembly the v3 definition is better
    global relative_empirical_first_zero_pos, relative_breakpoint_pos
    asm_penet = int(sat_end - sat_start +1)
    if sat_LR == "left":  # ref_point = sat_end 
        try:
            relative_empirical_first_zero_pos = str(int(ref_point - empirical_first_zero_pos))
        except:
            # print('ERROR in   relative_empirical_first_zero_pos = int(ref_point - empirical_first_zero_pos). relative_empirical_first_zero_pos value probably is None')
            relative_empirical_first_zero_pos = '.'        
        try:
            relative_x_intercept  = str(int(ref_point - x_intercept))
        except:
            print('ERROR in   relative_x_intercept  = int(ref_point - x_intercept). x_intercept value probably is None (artifact) ')
            # print('ref_point:',ref_point , '   x_intercept:',x_intercept)
            relative_x_intercept='.'
        relative_breakpoint_pos = int(closest_breakpoint - ref_point + 1)
    if sat_LR == "right": # ref_point = sat_start        
        try:
            relative_empirical_first_zero_pos = str(int(empirical_first_zero_pos - ref_point))
        except:
            # print('ERROR in   relative_empirical_first_zero_pos = str(int(empirical_first_zero_pos - ref_point)). relative_empirical_first_zero_pos value probably is None')
            relative_empirical_first_zero_pos = '.'           
        try:
            relative_x_intercept  = str(int(x_intercept - ref_point) )
        except:
            print('ERROR in   relative_x_intercept  = int(x_intercept - ref_point). x_intercept value probably is None (arifact)  ')
            # print('ref_point:',ref_point , '   x_intercept:',x_intercept)
            relative_x_intercept='.'            
        relative_breakpoint_pos = int(ref_point - closest_breakpoint + 1)                
    print(f"{sequenceID:<30}{monomer:>10}{sat_LR:>10}{ref_point:>10d}{asm_penet:>13d}{relative_empirical_first_zero_pos:>13s}{relative_x_intercept:>13s}{relative_breakpoint_pos:>10d}")    


    # Generate plot if requested
    if args.graph != 'none':
        # Save graph in current directory instead of data file directory
        # Get just the filename without path for the base name
        data_filename = os.path.basename(depth_data)
        base_name = os.path.splitext(data_filename)[0]
        
        # Build output filename with optional suffix
        if args.graph_suffix:
            # Remove any file extension from suffix if provided
            suffix = args.graph_suffix
            if suffix.startswith('.'):
                suffix = suffix[1:]
            output_file = f"{base_name}_fit{suffix}.{args.graph}"
        else:
            output_file = f"{base_name}_fit.{args.graph}"
        
        print(f"\nGenerating plot: {output_file}", file=log_file)
        
        # Pass ref_point and x_intercept info for plotting
        ref_point = ref_point if ref_point is not None else None
        create_plot(x_full_data, y_full_data, result, output_file, args.graph, 
                   sat_start, sat_end, ref_point, x_intercept, x_intercept_quadratic, 
                   breakpoint_ref_distance, empirical_first_zero_pos, args.graph_ref_lines)             


if __name__ == "__main__":
    main()

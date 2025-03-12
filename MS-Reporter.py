import configparser
import os
import re
import math
from collections import defaultdict

# -----------------------------------------------------------------------------------
# 1) READ CONFIGURATION FROM config.ini
# -----------------------------------------------------------------------------------

def read_config():
    config = configparser.ConfigParser()
    config_path = os.path.join(os.path.dirname(__file__), "config.ini")
    config.read(config_path)
    # Use fallback values if not defined in config.ini
    formula_str = config.get("DEFAULT", "formula")
    charge_range_default = config.get("DEFAULT", "preview_charge_range", fallback="1-30")
    output_mode = config.get("DEFAULT", "output_mode", fallback="decimals") # "decimals" or "sigfigs"
    precision = config.getint("DEFAULT", "precision", fallback=4)
    return formula_str, charge_range_default, output_mode, precision

# -----------------------------------------------------------------------------------
# 2) ISOTOPIC DATA (from NIST https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl)
# -----------------------------------------------------------------------------------
isotope_data = {
    'C': [(12.000000, 0.9893), (13.00335483507, 0.0107)],
    'H': [(1.00782503223, 0.999885), (2.01410177812, 0.000115)],
    'N': [(14.00307400443, 0.99636), (15.00010889888, 0.00364)],
    'O': [(15.99491461957, 0.99757), (16.99913175650, 0.00038), (17.99915961286, 0.00205)],
    'S': [(31.9720711744, 0.9499), (32.9714589098, 0.0075), (33.967867004, 0.0425), (35.96708071, 0.0001)],
    'B': [(10.01293695, 0.199), (11.00930536, 0.801)],
    'F': [(18.99840316273, 1)],
    'P': [(30.97376199842, 1)],
    'K': [(38.9637064864, 0.932581), (39.963998166, 0.000117), (40.9618252579, 0.067302)]
}

# Proton mass for ionization (in Da)
proton_mass = 1.007276466812

# -----------------------------------------------------------------------------------
# 3) HELPER FUNCTIONS FOR ISOTOPIC CALCULATIONS
# -----------------------------------------------------------------------------------

def parse_formula(formula_str):
    """
    Parse a molecular formula like "C383H619N121O110S3" 
    into a dictionary (e.g., {'C': 383, 'H': 619, ...})
    """
    pattern = r'([A-Z][a-z]*)(\d*)'
    formula = {}
    for (elem, count) in re.findall(pattern, formula_str):
        formula[elem] = formula.get(elem, 0) + (int(count) if count else 1)
    return formula

def merge_peaks(peaks, merge_threshold=0.001):
    """
    Merge peaks that are within merge_threshold (in Da) using weighted averaging.
    """
    if not peaks:
        return peaks
    peaks.sort(key=lambda x: x[0])
    merged = []
    current_mass, current_intensity = peaks[0]
    for mass, intensity in peaks[1:]:
        if mass - current_mass < merge_threshold:
            new_intensity = current_intensity + intensity
            current_mass = (current_mass * current_intensity + mass * intensity) / new_intensity
            current_intensity = new_intensity
        else:
            merged.append((current_mass, current_intensity))
            current_mass, current_intensity = mass, intensity
    merged.append((current_mass, current_intensity))
    return merged

def convolve(dist1, dist2, tol=1e-12, merge_tol=0.001):
    """
    Convolve two distributions (lists of (mass, intensity) pairs).
    """
    new_dict = {}
    for m1, p1 in dist1:
        for m2, p2 in dist2:
            prob = p1 * p2
            if prob < tol:
                continue
            new_mass = m1 + m2
            new_dict[new_mass] = new_dict.get(new_mass, 0) + prob
    new_peaks = list(new_dict.items())
    new_peaks = merge_peaks(new_peaks, merge_threshold=merge_tol)
    return new_peaks

def convolution_power(base, n, tol=1e-12, merge_tol=0.001):
    """
    Compute the n-fold convolution (power) of a base distribution using exponentiation by squaring.
    """
    if n == 0:
        return [(0.0, 1.0)]
    if n == 1:
        return base
    half = convolution_power(base, n // 2, tol, merge_tol)
    result = convolve(half, half, tol, merge_tol)
    if n % 2 == 1:
        result = convolve(result, base, tol, merge_tol)
    return result

def compute_isotopic_distribution(formula):
    """
    Given a formula dictionary, compute and return the isotopic distribution as a list of (mass, intensity).
    """
    distribution = [(0.0, 1.0)]
    for elem, count in formula.items():
        if elem not in isotope_data:
            raise ValueError(f"Isotopic data for element {elem} not provided.")
        elem_dist = convolution_power(isotope_data[elem], count)
        distribution = convolve(distribution, elem_dist)
    return distribution

def merge_for_instrument_resolution(distribution, fwhm=0.1):
    """
    Merge peaks using a threshold equal to FWHM (instrument resolution).
    """
    return merge_peaks(distribution, merge_threshold=fwhm)

def normalize_distribution(distribution):
    """
    Normalize intensities so that the most intense peak is 100%.
    """
    if not distribution:
        return distribution
    max_intensity = max(p for _, p in distribution)
    return [(mass, intensity / max_intensity * 100) for mass, intensity in distribution]

# -----------------------------------------------------------------------------------
# 4) EXPERIMENTAL PEAKS PARSING AND MATCHING
# -----------------------------------------------------------------------------------

def parse_experimental_peaks(input_text):
    """
    Parse experimental peak lines.
    Expected format: one peak per line with columns (index, m/z, Res., S/N, I %, FWHM)
    Returns a list of dictionaries with keys: 'mz' and 'I' (for intensity percentage).
    """
    peaks = []
    lines = input_text.strip().splitlines()
    for line in lines:
        # Skip header lines starting with '#' or if not enough columns.
        if line.strip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            mz = float(parts[1])
            intensity = float(parts[4])
            peaks.append({"mz": mz, "I": intensity})
        except:
            continue
    return peaks

def find_best_match(calc_mz, exp_peaks, output_mode, precision, ppm_accept=100, ppm_warn=20, ppm_search=1000):
    """
    For a given calculated m/z, find the best matching experimental peak.
    
    Returns: (selected_peak, list_of_notes)
    The calculated m/z is formatted using format_value().
    Experimental values in the notes are forced to 4 decimals.
    """
    notes = []
    candidates = []
    for peak in exp_peaks:
        ppm_diff = abs(peak["mz"] - calc_mz) / calc_mz * 1e6
        if ppm_diff <= ppm_search:
            candidates.append((peak, ppm_diff))
    if not candidates:
        notes.append(f"No experimental peak found within {ppm_search} ppm.")
        return None, notes
    # Filter candidates within the acceptance threshold (ppm_accept)
    candidates_accept = [(peak, ppm_diff) for (peak, ppm_diff) in candidates if ppm_diff <= ppm_accept]
    if not candidates_accept:
        notes.append(f"No experimental peak found within {ppm_accept} ppm for calculated m/z {format_value(calc_mz, output_mode, precision)}.")
        return None, notes
    # Select candidate by minimal ppm difference
    candidate_by_diff, diff_ppm = min(candidates_accept, key=lambda x: x[1])
    if diff_ppm > ppm_warn:
        notes.append(f"For calc m/z {format_value(calc_mz, output_mode, precision)}, the best match ({format_value(candidate_by_diff['mz'], output_mode, precision)}) is off by {diff_ppm:.1f} ppm (> {ppm_warn} ppm).")
    candidate_by_intensity, _ = max(candidates, key=lambda x: x[0]["I"])
    if abs(candidate_by_diff["mz"] - candidate_by_intensity["mz"]) > 1e-4:
        notes.append(f"Note: within 1000 ppm, the peak with highest I% is {candidate_by_intensity['mz']:.4f} (I%={candidate_by_intensity['I']}) which differs from the closest match {format_value(candidate_by_diff['mz'], output_mode, precision)}.")
    return candidate_by_diff, notes


# -----------------------------------------------------------------------------------
# 5) HELPER FUNCTIONS FOR FORMATTING OUTPUT
# -----------------------------------------------------------------------------------
def format_sigfig(num, sigfigs):
    """
    Format a number to the specified number of significant digits.
    """
    if num == 0:
        return "0." + "0" * (sigfigs - 1)
    order = math.floor(math.log10(abs(num)))
    decimals_needed = sigfigs - order - 1
    if decimals_needed < 0:
        decimals_needed = 0
    format_str = "{:." + str(decimals_needed) + "f}"
    return format_str.format(num)

def format_value(num, mode, precision):
    """
    Format a number according to the selected mode.
    mode: "decimals" or "sigfigs"
    If mode is "sigfigs", format the number to 'precision' significant figures.
    Otherwise, format the number with 'precision' fixed decimals.
    """
    if mode.lower() == "sigfigs":
        return format_sigfig(num, precision)
    else:
        return f"{num:.{precision}f}"

# -----------------------------------------------------------------------------------
# 6) MAIN CODE
# -----------------------------------------------------------------------------------

def main():
    # Read settings from config.ini
    formula_str, default_charge_range, output_mode, precision = read_config()
    
    print("Parsed formula from config.ini:", formula_str)
    print("Default protonation states range from config.ini:", default_charge_range)
    print("Output mode:", output_mode, "with precision", precision)
    print()
    
    # Parse the formula
    formula = parse_formula(formula_str)
    
    # Compute isotopic distribution
    print("Computing isotopic distribution (this may take a moment)...")
    raw_distribution = compute_isotopic_distribution(formula)
    merged_distribution = merge_for_instrument_resolution(raw_distribution, fwhm=0.1)
    norm_distribution = normalize_distribution(merged_distribution)
    norm_distribution.sort(key=lambda x: x[1], reverse=True)
    
    print("\nTop 3 isotopic masses (relative intensity):")
    for i, (mass, intensity) in enumerate(norm_distribution[:3], 1):
        print(f"{i}. Mass = {format_value(mass, output_mode, precision)} Da, Relative intensity = {intensity:.2f}%")
    
    most_abundant_mass = norm_distribution[0][0] if norm_distribution else 0.0

    # Show protonation states using the default range from config.ini.
    try:
        default_start, default_end = map(int, default_charge_range.split('-'))
    except Exception as e:
        print("Error reading default charge_range from config.ini. Using 1-30 as fallback.")
        default_start, default_end = 1, 30
    
    print(f"\nProtonation states (m/z) from +{default_start} to +{default_end}:")
    calc_states = {}
    for charge in range(default_start, default_end + 1):
        mz = (most_abundant_mass + charge * proton_mass) / charge
        calc_states[charge] = mz
        print(f"  [M+{charge}H] {charge}+  =>  {format_value(mz, output_mode, precision)} m/z")
    
    # Ask the user which protonation states to report.
    user_input = input(f"\nWhich protonation states do you want to report? (e.g. 7-10, press Enter for default {default_charge_range}): ").strip()
    if not user_input:
        user_input = default_charge_range

    try:
        start_str, end_str = user_input.split('-')
        start_charge = int(start_str)
        end_charge = int(end_str)
    except Exception as e:
        print("Invalid input. Use format '10-7' or '7-10'. Exiting.")
        return

    # Preserve the order as entered by the user
    if start_charge <= end_charge:
        charge_list = list(range(start_charge, end_charge + 1))
    else:
        charge_list = list(range(start_charge, end_charge - 1, -1))
    
    # Build final calculated report parts (plain ASCII)
    calc_report_parts = []
    for ch in charge_list:
        # If the requested charge state is not in calc_states, compute it.
        if ch not in calc_states:
            calc_states[ch] = (most_abundant_mass + ch * proton_mass) / ch
        mz = calc_states[ch]
        calc_report_parts.append(f"[M+{ch}H]{ch}+ " + format_value(mz, output_mode, precision))
    calc_report_parts.append("[M] " + format_value(most_abundant_mass, output_mode, precision))
    
    report_line = "m/z calcd for " + formula_str + " " + ", ".join(calc_report_parts)
    print("\nReport (plain ASCII):")
    print(report_line)
    
    # Ask if the user wants to compare with experimental masses
    compare_choice = input("\nDo you want to compare with experimental masses? (y/n): ").strip().lower()
    if compare_choice and compare_choice.startswith("y"):
        print("\nPlease paste the experimental peak list below (end input with an empty line):")
        exp_lines = []
        while True:
            line = input()
            if not line.strip():
                break
            exp_lines.append(line)
        exp_input = "\n".join(exp_lines)
        exp_peaks = parse_experimental_peaks(exp_input)
        if not exp_peaks:
            print("No valid experimental peaks were parsed.")
            return
        
        found_report_parts = []
        overall_notes = []
        comparison_parts = []
        # For each selected protonation state, find best experimental match.
        for ch in charge_list:
            calc_mz = calc_states[ch]
            best_match, notes = find_best_match(calc_mz, exp_peaks, output_mode, precision, ppm_accept=100, ppm_warn=20, ppm_search=1000)
            if best_match is None:
                found_report_parts.append("N/A")
                comparison_parts.append(f"[M+{ch}H]{ch}+ calc {calc_mz:.4f}, exp N/A")
            else:
                found_report_parts.append(format_value(best_match['mz'], output_mode, precision))
                error = abs(calc_mz - best_match['mz']) / calc_mz * 1e6
                # Force 4 decimals in the comparison output regardless of config value.
                comparison_parts.append(f"[M+{ch}H]{ch}+ calc {calc_mz:.4f}, exp {best_match['mz']:.4f} ({error:.1f} ppm)")
            if notes:
                overall_notes.append(f"For [M+{ch}H]{ch}+: " + " ".join(notes))
        
        # Also process the neutral (intact mass) [M]
        neutral_calc = most_abundant_mass
        neutral_match, neutral_notes = find_best_match(neutral_calc, exp_peaks, output_mode, precision, ppm_accept=100, ppm_warn=20, ppm_search=1000)
        if neutral_match is None:
            found_neutral = "N/A"
            comp_neutral = f"[M] calc {neutral_calc:.4f}, exp N/A"
        else:
            found_neutral = format_value(neutral_match['mz'], output_mode, precision)
            neutral_error = abs(neutral_calc - neutral_match['mz']) / neutral_calc * 1e6
            comp_neutral = f"[M] calc {neutral_calc:.4f}, exp {neutral_match['mz']:.4f} ({neutral_error:.1f} ppm)"
        found_report_parts.append(found_neutral)
        comparison_parts.append(comp_neutral)
        if neutral_notes:
            overall_notes.append("For [M]: " + " ".join(neutral_notes))

        # Build full report line: calculated masses and found masses
        full_report_line = "m/z calcd for " + formula_str + " " + ", ".join(calc_report_parts) \
                           + ", found m/z " + ", ".join(found_report_parts)
        print("\n" + full_report_line)
        print()  # Empty line between found and comparison lines
        comparison_line = "Comparison: " + ", ".join(comparison_parts)
        print(comparison_line)
        if overall_notes:
            print("\nNotes:")
            for note in overall_notes:
                print(" - " + note)
    else:
        print("\nNo experimental comparison requested.")

if __name__ == "__main__":
    main()

MS-Reporter
===========

Description:
------------
MS-Reporter is a Python program that calculates the isotopic distribution and theoretical m/z 
values for a given molecule based on its molecular formula. It also allows you to compare 
the calculated m/z values with experimental mass spectrometry data.

Features:
---------
- Computes the isotopic distribution using standard NIST isotopic data.
- Calculates theoretical m/z values for the intact molecule ([M]) and for a range of protonation 
  states (e.g., [M+H]+, [M+2H]2+, etc.).
- Reads configuration settings from a config.ini file.
- Optionally compares the calculated m/z values with experimental data, reporting the best match 
  and the ppm error for each protonation state (and for [M]).
- Reports all values in plain ASCII text.

Configuration:
--------------
A file named "config.ini" should be placed in the same folder as the Python script.
Below is an example of the config.ini file:

```
[DEFAULT]
formula = C675H1098N198O203S6
preview_charge_range = 5-30
decimals = 4
```

- "formula"      : The molecular formula (e.g., C675H1098N198O203S6).
- "preview_charge_range": The default range of protonation states to preview, in the format start-end (e.g., 5-30).
- "decimals"     : The number of decimal places for the calculated m/z values (except for the comparison, which always uses 4 decimals).

Usage:
------
1. Place MS-Reporter.py and config.ini in the same folder.
2. Edit config.ini as needed.
3. Run the program from a terminal (or command prompt):
      python MS-Reporter.py

Program Workflow:
-----------------
1. The program reads the molecular formula, the default protonation state range, and the 
   number of decimals from config.ini.
2. It computes the isotopic distribution for the provided formula and prints the top 3 
   isotopic masses along with their relative intensities.
3. It displays the calculated m/z values for all protonation states within the default range.
4. You are prompted to enter the protonation states you want to report. If you press Enter, 
   the default range from config.ini is used.
5. A report line is printed showing the calculated m/z values for the selected protonation states 
   and the intact molecule ([M]).
6. You are then asked if you wish to compare the calculated m/z values with experimental data.
   - If you answer "y", you will be prompted to paste your experimental peak list.
   - The program expects each line of the experimental peak list to contain columns: 
     index, m/z, Res., S/N, I %, FWHM (lines starting with '#' are ignored).
   - For each selected protonation state (and for [M]), the program finds the best matching 
     experimental peak within defined ppm thresholds and prints:
        a) A "found m/z" line showing the experimental values.
        b) A "Comparison:" line that lists for each state the calculated m/z (4 decimals), 
           the experimental m/z (4 decimals), and the ppm error.
        c) Any notes if the match is not ideal.
7. If you choose not to compare, only the calculated values are reported.

Example Output:
---------------
```Parsed formula from config.ini: C675H1098N198O203S6
Default protonation states range from config.ini: 5-30
Using 4 decimals for the output.

Computing isotopic distribution (this may take a moment)...
Top 3 isotopic masses (relative intensity):
1. Mass = 15427.0000 Da, Relative intensity = 100.00%
2. Mass = 15426.1000 Da, Relative intensity = 93.45%
3. Mass = 15425.2000 Da, Relative intensity = 85.67%

Protonation states (m/z) from +5 to +30:
  [M+5H] 5+  =>  3085.4000 m/z
  [M+6H] 6+  =>  2571.2000 m/z
  ...
  [M] 15427.0000 m/z

Which protonation states do you want to report? (e.g. 7-10, press Enter for default 5-30): 7-10

Report (plain ASCII):
m/z calcd for C675H1098N198O203S6 [M+7H]7+ 2204.9000, [M+8H]8+ 1929.4000, [M+9H]9+ 1715.1000, [M+10H]10+ 1543.7000, [M] 15427.0000

Do you want to compare with experimental masses? (y/n): y

Please paste the experimental peak list below (end input with an empty line):
#   m/z     Res.    S/N     I %     FWHM
1   2204.9  25000   15.0    35.0    0.0
2   1929.4  24000   14.0    32.0    0.0
3   1715.1  23500   13.5    30.0    0.0
4   1543.7  23000   12.5    28.0    0.0
5   15427.0 22000   20.0    50.0    0.0

found m/z 2204.9000, 1929.4000, 1715.1000, 1543.7000, 15427.0000
Comparison: [M+7H]7+ calc 2204.9000, exp 2204.9000 (0.0 ppm), [M+8H]8+ calc 1929.4000, exp 1929.4000 (0.0 ppm), [M+9H]9+ calc 1715.1000, exp 1715.1000 (0.0 ppm), [M+10H]10+ calc 1543.7000, exp 1543.7000 (0.0 ppm), [M] calc 15427.0000, exp 15427.0000 (0.0 ppm)

Notes:
 - For [M+9H]9+: For calc m/z 1715.1000, the best match (1715.1000) is off by 25.3 ppm (> 20 ppm).
```

Troubleshooting:
----------------
- Ensure your experimental peak list has at least five columns per line (index, m/z, Res., S/N, I %, FWHM).
- Lines starting with '#' are ignored.
- PPM thresholds for matching (default: 100 ppm acceptance, 20 ppm warning, 1000 ppm search) can be adjusted in the source code.


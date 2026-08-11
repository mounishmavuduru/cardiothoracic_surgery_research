"""Contextual check: a numeral can be present in the master and still be
attached to the wrong subject after paraphrase. For each substantive claim in
the condensed article, require the number to appear in the master within a
window that also contains the claim's subject keyword."""
import re
import sys

master = open(sys.argv[1], encoding="utf-8").read()
flat = " ".join(master.split()).replace("{,}", "").replace("\\,", "")

# (number-as-written, subject regex that must co-occur within +/-WINDOW chars)
WINDOW = 420
CLAIMS = [
    ("0.012", r"combined|gbt|gradient"),
    ("-0.007", r"combined|logistic|lr"),
    ("0.884", r"combined|baseline|competitor"),
    ("0.851", r"combined|gradient|gbt"),
    ("0.449", r"0.012|combined"),
    ("0.70", r"-0.007|combined"),
    ("-0.021", r"Roney"),
    ("0.009", r"Roney"),
    ("0.035", r"UW|Boyle"),
    ("0.870", r"Roney"),
    ("0.834", r"UW|Boyle|logistic"),
    ("0.748", r"UW|Boyle|gradient"),
    ("34", r"Roney|inducible"),
    ("182", r"combined|real|meshes"),
    ("42", r"inducible|combined"),
    ("0.0028", r"subspace"),
    ("4150", r"power|cases|subspace"),
    ("23000", r"single-vector|power|cases"),
    ("0.0036", r"Kuramoto"),
    ("0.38", r"Kuramoto"),
    ("0.12", r"FitzHugh|FHN"),
    ("0.996", r"canonical"),
    ("0.289", r"permut"),
    ("0.970", r"R\^\{?2|reconstruct|regress"),
    ("2422", r"rho|\\rhoval|median|validity"),
    ("1082", r"interquartile|IQR|\["),
    ("16554", r"range"),
    ("272", r"factor|smallest"),
    ("1.00", r"Weyl|exponent|N\^"),
    ("0.5\\%", r"Weyl|exponent"),
    ("0.117", r"grouping|naive|shape-family"),
    ("0.640", r"grouped|shape-family"),
    ("0.756", r"naive|random"),
    ("55000", r"network|p-value|DeLong"),
    ("0.197", r"clinical|logistic|SFI"),
    ("0.250", r"baseline|anti-predictive"),
    ("0.447", r"SFI|below chance"),
    ("0.0061", r"DeLong|clinical"),
    ("0.137", r"noise"),
    ("0.333", r"empirical|noise"),
    ("0.995", r"prevalence|predicted probability"),
    ("0.345", r"within folds|baseline"),
    ("48", r"recur|event"),
    ("0.24", r"power"),
    ("371", r"patients|power"),
    ("859", r"patients|power"),
    ("0.435", r"univariate|chance"),
    ("0.547", r"univariate|chance"),
    ("58", r"1500|resolution|inducibility rate"),
    ("37.5", r"50|resolution|inducibility rate"),
    ("34.0", r"2000|operating|Roney"),
    ("7/24", r"consistent|tier"),
    ("650", r"reentry|self-sustained"),
    ("150", r"burst|cycle length"),
    ("0.36", r"Delta w|uncoupling"),
    ("0.3", r"c_|perp|1\\%"),
    ("0.05", r"varepsilon|gate|threshold|learning rate"),
    ("120", r"trees"),
    ("0.8", r"subsample|source study"),
    ("0.603", r"fibrosis-only|degrad"),   # results-only; expect MISS in master
    ("0.537", r"gradient|gbt"),           # results-only; expect MISS in master
]

miss = []
for num, subj in CLAIMS:
    pat = re.escape(num.replace("\\%", "%")).replace(r"\-", "-")
    ok = False
    for m in re.finditer(r"(?<![\d.])" + pat + r"(?![\d])", flat):
        ctx = flat[max(0, m.start() - WINDOW):m.end() + WINDOW]
        if re.search(subj, ctx, re.I):
            ok = True
            break
    if not ok:
        miss.append((num, subj))

print(f"claims checked: {len(CLAIMS)}   unsupported in master: {len(miss)}")
for num, subj in miss:
    print(f"  {num}  (no master context matching /{subj}/)")

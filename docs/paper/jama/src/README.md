# Sources for the condensed JAMA Cardiology Original Investigation

`JAMACardio_OriginalInvestigation.tex` is generated, not hand-edited. Edit the
fragments here and regenerate, or the next build silently discards the change.

## Build

```
python docs/paper/jama/src/assemble_oi.py docs/paper/jama/JAMACardio_OriginalInvestigation.tex
tools/tectonic.exe -X compile docs/paper/jama/JAMACardio_OriginalInvestigation.tex \
    --outdir docs/paper/jama --keep-logs
```

`assemble_oi.py` prints the word counts and substitutes them into the title
page, so the declared counts cannot drift from the text. It aborts without
writing if a citation-injection pattern fails to match exactly once.

## Verification

```
python docs/paper/jama/src/verify_oi.py \
    docs/paper/jama/JAMACardio_OriginalInvestigation.tex docs/paper/manuscript.tex
python docs/paper/jama/src/claims_oi.py docs/paper/manuscript.tex
```

`verify_oi.py` requires every numeral in the article body to appear either in
the audited master manuscript or in a `results/*.json` file at some legitimate
rounding. `claims_oi.py` is the stronger check: it requires each substantive
number to occur in the master *near its own subject*, which catches a numeral
that survived paraphrase but changed what it refers to.

Known and intentional `claims_oi.py` reports:

* `0.603`, `0.537` --- clinical gradient-boosting and fibrosis-only figures that
  the master states in prose without the numeral. Both are traced to
  `results/clinical_endpoint.json`
  (`baseline_degradation_section_8_6.fibrosis_only_auc.gbt` = 0.6029,
  `primary_section_8_4.result.gbt.grouped_auc_sfi` = 0.5368).
* `0.5\%` --- present in the master as a LaTeX-escaped percent sign, which the
  checker's literal pattern does not match (master line 1259).

## Limits enforced

| item | JAMA limit | this article |
|---|---|---|
| main text | 3000 words | 2980 |
| abstract | 350 words, structured | 350, seven headings |
| Key Points | 75--100 words | 95 |
| tables + figures | 5 | 5 (2 tables, 3 figures) |
| references | 50--75 typical | 19 |

The reference count is below the range JAMA describes as typical. It is not
padded: the argument cites the identities, the datasets, the statistical
methods, and the inducibility literature it actually uses. Supplement 2 carries
the full technical report.

## Relationship to the other files here

* `docs/paper/manuscript.tex` --- the audited master, source of truth for every
  number.
* `Supplement2_full_technical_report.tex` --- the master with its title block
  replaced, built by `build_supp2.py`. Body byte-identical to the master.
* `Supplement1_TRIPOD_AI_checklist.tex` --- 52 sub-items, scoped to the clinical
  endpoint only.

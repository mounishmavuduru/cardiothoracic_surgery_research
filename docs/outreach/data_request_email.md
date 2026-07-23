# Data-request email — CEMRG / Roney group (per-patient recurrence labels)

*Fill the [BRACKETS] before sending. Primary target: the Roney/CEMRG group, who release 100 public
LA models and demonstrably hold the linked 1-year recurrence outcomes (34/99 recurred). A UW/Boyle
variant is below. Keep it short, honest, specific, and low-effort for them to say yes to.*

**To:** [corresponding author — e.g., Dr Caroline Roney / the CEMRG contact listed on cemrg.com]
**Cc:** [your mentor / supervising teacher, if any]
**Subject:** Student research request: per-patient recurrence labels for the 100 public LA models (CIRCEP 2022)

---

Dear Dr [Last name],

I'm [Full name], a [grade / year] student at [School], carrying out an independent computational
research project on atrial-fibrillation substrate analysis (supervised by [mentor, if any]). I've been
working with your publicly released 100 patient-specific left-atrial models (CEMRG; *Circ Arrhythm
Electrophysiol* 2022, doi:10.1161/CIRCEP.121.010253), which are an excellent, analysis-ready resource —
thank you for making them open.

My project evaluates whether spectral graph-theoretic features of the atrial conduction network
(algebraic connectivity and its edge-sensitivities) add predictive value beyond established fibrosis
and connectivity measures. So far I have tested this only against an in-silico inducibility labeller,
which I know is a limitation. To make the evaluation meaningful, I would like to test against a **real
clinical outcome**.

Would it be possible to share the **per-patient 1-year post-ablation recurrence labels** (the outcome
used in the 2022 study, keyed to the 100 released model IDs)? A single de-identified column
(model ID → recurred yes/no, and blanking-period/censoring if available) is all I would need; I do not
need any additional patient information.

In return I'm glad to: (i) share my full analysis code and results with your group; (ii) use the data
only for this non-commercial student project under whatever data-use terms you specify; and (iii)
acknowledge your group (or include you as appropriate) in any resulting write-up. I'm also happy to
sign a data-use agreement or route the request through your institution's process if that's cleaner.

I completely understand if the outcome labels can't be shared for confidentiality reasons — in that
case, any pointer to a public dataset that links atrial geometry to real outcomes would be hugely
appreciated.

Thank you very much for your time and for the open resource.

Best regards,
[Full name]
[School] · [City, Country]
[email] · [optional: project one-pager / GitHub link]

---

## UW/Boyle variant (secondary)

Same body, retargeted:
- **To:** [Prof. Patrick Boyle / corresponding author, University of Washington]
- **Subject:** Student request: per-patient recurrence labels for the Dryad atrial-mesh cohort (Comms Med 2025)
- Swap the dataset sentence for: *"I've been working with your publicly released atrial meshes
  (Dryad doi:10.5061/dryad.kkwh70sg0; Communications Medicine 2025, doi:10.1038/s43856-025-01058-4),
  82 patients pre/post ablation."* and request the per-patient recurrence labels keyed to IDs 001–015
  and 021–087.

## Honest expectations (read before sending)
- This is a **lottery ticket**, not a plan. Groups often can't share outcome labels even when willing,
  and replies can take weeks. Send it, then proceed with the rigorous (Path B) manuscript regardless.
- If labels *do* arrive: the bar for a novel result is **beating their published predictor** (their RF
  reached AUROC ~0.80) using spectral features — a real but hard target that may still come up null.

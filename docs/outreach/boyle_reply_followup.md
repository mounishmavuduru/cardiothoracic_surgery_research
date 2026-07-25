# Reply to Prof. Boyle — thanks + four questions that gate the analysis

*Drafted 2026-07-24, after the recurrence outcomes arrived. Fill the [BRACKETS] before
sending. Keep it short: he has already done us a large favour, and each question below is
answerable in one line.*

**To:** Prof. Patrick M. Boyle <pmjboyle@uw.edu>
**Subject:** Re: recurrence labels — thank you, and four quick clarifications

---

Dear Prof. Boyle,

Thank you — genuinely. The outcome column is the single thing my project could not obtain
any other way, and I know sharing it took effort on your end.

Before I touch it, I have frozen and timestamped the analysis plan so the test stays
falsifiable: endpoints, feature sets, multiplicity correction, and a power calculation done
from the marginal event counts alone. It is public at [REPO URL], git-tagged
`prereg-real-outcomes-20260724`. Two things I want to flag honestly up front:

- With 82 patients and 48 events, the minimum effect I can detect is roughly ΔAUC 0.10–0.15.
  I have pre-declared that a smaller null will be reported as *inconclusive*, not as evidence
  against the method.
- I am **not** attempting to beat your 0.80 AUC model. You have 89 features including EHR
  risk factors; I have mesh geometry only. Any comparison I publish will say so plainly.

Four clarifications would materially improve the work:

1. **Redistribution.** My default is that the outcome column is restricted-use: it stays out
   of my public repository, and only aggregate statistics appear in any write-up. Is that the
   right reading of your terms, or is there a form of sharing you would prefer? I am happy to
   sign a DUA or route this through UW's process.

2. **Confirming the `elemTag` values.** The README does not define them, so I worked them out
   geometrically. Element fractions across all 164 meshes:

   | tag | pre-ablation | post-ablation | my reading |
   |---|---|---|---|
   | 111 | 59.0 % | 47.3 % | healthy myocardium |
   | 115 | 20.9 % | 14.9 % | fibrosis |
   | 164 | 20.0 % | 19.9 % | the caps over the four PVs and the mitral valve |
   | 199 | absent | 17.9 % | ablation lesion |

   I am fairly confident about 164 specifically: it forms exactly five connected components,
   each topologically a disc, 11–42 mm across; it has *zero* edge-adjacency to tag 115
   anywhere; and deleting it leaves a surface of Euler characteristic −3, i.e. exactly five
   boundary loops. Could you confirm? I ask because my loader had been treating 164 as
   half-fibrotic tissue, which not only added ~20 % spurious fibrosis to every patient but
   sealed the atrium's orifices, so activation could cross the mitral valve rather than
   circle it.

3. **Continuous fibrosis.** Is the underlying continuous LGE intensity-ratio (IIR) field per
   vertex something you could share? A binary fibrotic/not tag weakens the fibrosis baseline
   I am trying to beat, which biases my comparison in my own favour — I would rather remove
   that confound than disclose it.

4. **The released fibre field is degenerate — this one may matter to you.** Each mesh does
   contain a `VECTORS fiber` array, but in all 164 files it is a constant `(1, 0, 0)` for
   every element: mean directional spread is exactly 0. I suspect fibres did not survive the
   downsampling to 0.5 mm average edge length described in the README. It matters because
   anyone following the README's suggestion to run reaction–diffusion simulations on these
   meshes gets globally uniform anisotropy rather than atrial fibre architecture. Would it be
   possible to share the fibre field, or is there a re-export that preserves it?

5. **Follow-up completeness.** For the `NR` patients, were all of them followed the full two
   years, or are some censored early / lost to follow-up? And was a standard 90-day blanking
   period applied before counting recurrence? This changes whether I treat the outcome as
   binary or time-to-event.

Whatever the result — and a null is a live possibility given the power above — I will share
the code and findings with you before anything goes public, and I would be glad to
acknowledge or include you and your group as you see fit.

Thank you again for the meshes, the paper, and the labels.

Best regards,
[Full name]
[School] · [City]
[email] · [repo link]

---

## Why each question matters (internal notes, do not send)

- **Q1** is the one that can cause real harm if I guess wrong. `docs/DATASETS.md` recorded
  these labels as withheld "protecting confidential patient information", so the CC0 covering
  the meshes clearly does not extend to them. Default to non-redistribution until he says
  otherwise, in writing.
- **Q2 is the real blocker, and it is a likely bug in our own loader.**
  `src/asb/substrate/uw_boyle.py:68` maps `elemTag` → fibrosis as
  `{111:0.0, 115:1.0, 164:0.5, 199:1.0}`. The measured tag fractions say 164 occupies
  20.0 % of cells pre-ablation and 19.9 % post — i.e. **ablation does not touch it**, which is
  not how fibrosis behaves and is exactly how a valve annulus or PV sleeve behaves. If that is
  right, the 0.5 mapping adds a near-constant ~20 % pseudo-fibrosis to every patient. That
  does two bad things: it inflates `fibrosis_burden` while *shrinking its between-patient
  variance* (weakening the competitor baseline), and it corrupts the Δw uncoupling field that
  SFI itself is built from. Do not run the confirmatory analysis until this is resolved.
- **Q3** removes the last substrate-fidelity confound. Note the cohort-split question was
  dropped: the Dryad README already confirms `ID001–015` = holdout, `ID021–087` = original
  cohort, so secondary analysis §8.5.1 can honour their split without asking.
- **Q4** decides binary logistic vs. Cox. If NR includes patients lost at, say, 8 months, the
  binary outcome is biased toward NR and every AUC is optimistic in an uncontrolled way.

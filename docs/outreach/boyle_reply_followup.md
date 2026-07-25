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

2. **Fibrosis fidelity.** The Dryad `.vtk` files carry a per-element `elemTag` (I read
   111 / 115 / 164 / 199 as healthy / fibrotic / boundary-or-scar). My method needs a
   *continuous* fibrosis field, and using the categorical tag measurably weakens the fibrosis
   baseline I am trying to beat — which would bias my result in my own favour. Is the
   underlying continuous LGE intensity-ratio (IIR) field per vertex something you could share,
   or is the tag the intended fidelity for these released meshes?

3. **Cohort split.** Am I right that `ID001–ID015` is the holdout and `ID021–ID087` the
   training cohort from the paper, and that IDs 16–20 were excluded upstream? I would like to
   honour your original split rather than invent my own.

4. **Follow-up completeness.** For the `NR` patients, were all of them followed the full two
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
- **Q2** is the technical blocker. Per `src/asb/substrate/uw_boyle.py:68` the loader maps
  `elemTag` to fibrosis in {0.0, 0.5, 1.0}. Compared with the Roney cohort's continuous `IIR`
  ramp, that is a coarse substrate. It degrades the competitor baseline more than it degrades
  SFI, which inflates ΔAUC in our favour — pre-registration §8.6 discloses this and commits to
  reporting it as a confound. A continuous IIR field would remove the confound entirely.
- **Q3** costs nothing to ask and lets secondary analysis 8.5.1 be a genuinely external test.
- **Q4** decides binary logistic vs. Cox. If NR includes patients lost at, say, 8 months, the
  binary outcome is biased toward NR and every AUC is optimistic in an uncontrolled way.

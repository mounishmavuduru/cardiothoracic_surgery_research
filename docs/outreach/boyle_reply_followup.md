# Reply to Prof. Boyle — ready to send

*Rewritten 2026-07-25. Every factual claim below has been checked against a measurement made
this session; see the verification notes at the bottom. Fill `[REPO URL]` before sending.*

**To:** Prof. Patrick M. Boyle <pmjboyle@uw.edu>
**Subject:** Re: recurrence labels — thank you, and a few questions (plus one issue with the Dryad meshes)

---

Dear Prof. Boyle,

Thank you — sincerely. The outcome column was the one thing I had no other way to obtain, and
I know that sharing it took effort on your end.

Before touching it I froze and timestamped the analysis plan, so the test stays falsifiable:
endpoints, feature sets, multiplicity correction, and a power calculation done from the
marginal event counts alone. It is public at [REPO URL], git-tagged
`prereg-real-outcomes-20260724`. Two things I want to be upfront about:

- With 82 patients and 48 events, the smallest effect I can detect is roughly ΔAUC 0.10–0.15.
  At ΔAUC 0.05 I have about 25 % power. I have pre-declared that anything smaller will be
  reported as inconclusive rather than as evidence against the method.
- I am not trying to beat your model. It draws on EHR risk factors and simulation-derived
  features that I do not have; mine uses mesh geometry alone. Any comparison I publish will
  say so plainly.

Five questions, each answerable in a line or two:

**1. Redistribution.** My default is that the outcome column is restricted-use: it stays out
of my public repository, and only aggregate statistics appear in any write-up. Is that the
right reading of your terms, or would you prefer something different? I am happy to sign a
DUA or route this through UW's process.

**2. `elemTag` semantics.** The README does not define the values, so I worked them out
geometrically. Measured across all 164 meshes:

| tag | pre-ablation | post-ablation | my reading |
|---|---|---|---|
| 111 | 59.0 % | 47.3 % | healthy myocardium |
| 115 | 20.9 % | 14.9 % | fibrosis |
| 164 | 20.0 % | 19.9 % | caps over the pulmonary veins and mitral valve |
| 199 | absent | 17.9 % | ablation lesion |

I am most confident about 164. It resolves into four to six large components per mesh — five
in 75 of the 82, four in five of them and six in two, which is about what I would expect from
normal pulmonary-vein variation — each component is topologically a disc, and it borders
fibrotic elements on only ~0.08 % of its edges against a 14–20 % chance rate. Removing it
opens the surface. Have I read that correctly? I ask because I had been treating 164 as
half-fibrotic tissue, which sealed the atrium's orifices and let activation cross the mitral
valve instead of circling it.

**3. Continuous fibrosis.** Would it be possible to share the per-vertex LGE intensity ratio?
A binary fibrotic/not tag weakens the fibrosis baseline that my method has to beat, which
biases the comparison in my own favour — I would rather remove that confound than disclose it.

**4. The fibre field looks as though it did not survive export.** Every mesh contains a
`VECTORS fiber` array, but in all 164 files it is a constant `(1, 0, 0)`: the mean directional
spread is exactly zero. I suspect it was lost in the downsampling to 0.5 mm average edge
length. I mention it because the README suggests these meshes for reaction–diffusion work, and
anyone doing that would get globally uniform anisotropy rather than atrial fibre architecture.
For what it is worth, I checked whether it actually mattered for my own results by stripping
the fibres from a different cohort that has real ones, and the inducibility rate did not move
at all — so this is not me attributing my results to your data. Is there a re-export that
preserves the fibres?

**5. Follow-up completeness.** Were all of the `NR` patients followed the full two years, or
are some censored early? And was a blanking period applied before recurrence was counted?
This decides whether I treat the outcome as binary or as time-to-event.

Whatever the result — and a null is a live possibility given the power above — I will share
the code and findings with you before anything goes public, and I would be glad to acknowledge
or include you and your group as you see fit.

Thank you again for the meshes, the paper, and the labels.

Best regards,
Mounish Mavuduru

---

## Verification notes (internal — do not send)

Every number in the email traces to a measurement made this session:

- **Tag fractions** — census over all 164 meshes.
- **"four to six large components; five in 75 of 82"** — connected components of the tag-164
  region at ≥1 % of the region, all 82 pre-ablation meshes: 5 in 75, 4 in 5, 6 in 2.
- **"topologically a disc"** — Euler characteristic χ = 1 per large component.
- **"~0.08 % of its edges against 14–20 % by chance"** — median 27 contacting edges; chance
  level is the fibrotic fraction.
- **"removing it opens the surface"** — χ falls by a median of 7 (range 2–12). Deliberately
  *not* phrased as "exactly five openings": an earlier draft said that, and it was wrong,
  generalised from mesh ID001 alone.
- **"constant (1,0,0), mean directional spread exactly zero"** — all 164 meshes.
- **"stripping the fibres… the rate did not move at all"** — Roney control, 20/62 with real
  fibres and 20/62 with constant fibres, McNemar p = 1.0.
- **Power figures** — `scripts/power_real_outcomes.py`, from the marginal 48/34 split only.

Claims deliberately **removed** from the earlier draft because measurement refuted them:

- "exactly five connected components" → it is 4–6.
- "zero edge-adjacency to fibrosis" → it is ~0.08 %, not zero.
- "Euler characteristic −3, i.e. exactly five boundary loops" → true of ID001 only.
- Any suggestion that the degenerate fibre field explains our low inducibility → the Roney
  control refutes it, and the email now says so explicitly.

The cohort-split question was dropped: the Dryad README already states that ID001–015 are the
holdout patients and ID021–087 the original cohort, so asking would have looked careless.

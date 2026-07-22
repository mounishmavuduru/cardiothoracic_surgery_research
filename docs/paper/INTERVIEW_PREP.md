# Interview preparation — hardest questions, honest answers

The questions most likely to sink the project, and the answer that survives an expert. Every answer
is defensible from `results/*.json`. The rule: never overclaim, always give the number.

---

**Q. What is *your* contribution? The math is all classical.**
Correct — and I say so first, before you do. Fiedler sensitivity (Ghosh–Boyd 2006), the master
stability function (Pecora–Carroll 1998), Weyl's law, Davis–Kahan — all classical. My contribution is
not mathematics. It is: (1) a pre-registered, leakage-controlled falsification of a plausible
biomarker on two real cohorts and 100,000 networks; (2) using ρ = ‖ΔL‖/(λ₃−λ₂) as an *a priori
validity screen* for spectral biomarkers and measuring that real tissue violates it by 10³; (3) a
Weyl scaling law showing the failure is lawful, not incidental; (4) a demonstration that standard
methodology would have called this a hit. The novelty is in the *use and the rigor*, not the theorems.

**Q. So your biomarker didn't work. Why should I care about a negative result?**
Two reasons. First, I pre-registered the null as an accepted outcome *before seeing any labels*, so
this is a clean falsification, not a fishing expedition — and the field is full of spectral biomarkers
that were never tested this way. Second, I don't just report "it failed"; I explain *why any biomarker
of this class must fail on real excitable media* (the validity radius) and give a one-line test anyone
can run before trusting such a marker. A negative result with a mechanism and a general precondition
is more useful than another over-claimed positive.

**Q. You got p < 10⁻⁴ for the subspace SFI at large N. So it works, right?**
No — and this is the most important number in the project. p < 10⁻⁴ only because I have 55,000+
networks. The *effect size* is ΔAUC = +0.003 on a base of 0.91. To detect that in humans at 80% power
needs ~4,150 patients; the largest real AF imaging cohorts are ~1,000. It is statistically significant
and clinically meaningless. That gap between significance and effect size is exactly the trap this
project is built to expose: with enough data, useless effects look real. Effect size, not p, decides.

**Q. What is your single weakest point?**
The labels are a monodomain simulator's inducibility verdict, not clinical outcomes and not the
openCARP gold standard. I'm explicit about it everywhere. It's a limitation of fidelity, not of logic:
the *methodology* (pre-registration, grouping, effect-size gate, spatial null, cross-medium transfer)
is label-agnostic, and the validity-radius argument is analytical and independent of the labeller.
The deferred openCARP pass would raise label fidelity, not change the ρ ≈ 2422 verdict.

**Q. Why should synthetic networks tell me anything about a real heart?**
They aren't a substitute for the heart — the real cohorts (Roney n=62, UW/Boyle n=82) are. The
synthetic networks do two things the real cohorts can't: they let me push N to 100,000 so the effect
size and its confidence interval converge (which is how I can state "+0.003, needs 4,150 cases"), and
they let me test whether the *pattern transfers to other excitable and oscillator media* — which
isolates what's a property of the mathematics versus an accident of one dataset. The real cohorts
carry the anatomy; the synthetics carry the scaling and the generality.

**Q. Fibrosis also localizes the origin. So what does |∇φ₂| add?**
For pure localization on real tissue, honestly, not much extra — and I show *why*: in a diffusively
coupled medium the Fiedler bottleneck *is* the low-conduction lesion, so the connectivity localizer
and the substrate localizer point at the same place by mechanism, not coincidence. What |∇φ₂| adds is
(a) it is purely connectivity-derived and still beats the graph-centrality baselines (Perron, weighted
degree), and (b) in the oscillator medium it localizes more sharply than the raw substrate. I tried to
find a regime where it beats imaging outright (a healthy source–sink isthmus); whether that works is an
open experiment I report either way.

**Q. Isn't the validity-radius section just "perturbation theory breaks when the perturbation is big"?**
Yes — and I state exactly that in the paper. It is textbook. I do not call it a theorem. What's mine is
turning it into an operational screening test for biomarkers (compute ρ; if ρ ≫ 1, distrust the linear
spectral marker), measuring ρ ≈ 2422 on real atria, and showing via a Weyl scaling law that ρ grows
lawfully with mesh resolution so the failure is structural. If you thought I was claiming a new theorem,
I've mis-written it; the claim is a framing-and-measurement contribution.

**Q. Did you rig the dissociation experiment to make the connectivity localizer win?**
No — and it *didn't* win, which is the tell. I built a knob (fibrosis-to-bottleneck distance) and let
the simulator decide the origin. In the abstract graph models it doesn't dissociate: when I move
fibrosis off the bottleneck the instability stops happening, because the bottleneck only becomes
vulnerable when low conduction makes it so — and that's what imaging sees. I report that negative and
the mechanism it revealed. The one model that could dissociate them (wavefront-curvature source–sink on
a healthy isthmus) is a harder, honest test I'm running separately.

**Q. How would this generalize beyond cardiology?**
The validity screen is domain-general. Spectral-graph fragility markers are used for seizure-focus
localization, brain connectomics, and power-grid vulnerability. The same precondition applies: compute
ρ = ‖ΔL‖/(gap); if it's ≫ 1, the first-order spectral marker is outside its validity radius and its
in-sample correlations can't be trusted. That's a portable, testable claim — arguably more useful than
a single-domain positive would have been.

**Q. What would you do next, with more time and compute?**
Three things: (1) the openCARP gold-standard labelling to raise fidelity (deferred for compute);
(2) finish the healthy-isthmus source–sink experiment to settle whether connectivity ever beats imaging;
(3) test the ρ screen on a genuinely different real domain (e.g., iEEG seizure networks) to demonstrate
the precondition transfers to a non-cardiac biomarker.

---

**Meta-note for the presenter.** Your strongest move in the room is to *volunteer* the limitations
before the judge finds them, and to lead with the validity radius and the p-value/effect-size
distinction rather than the biomarker. You are not selling a discovery; you are demonstrating that you
can kill your own hypothesis correctly and explain why. That is what graduate-level research looks like.

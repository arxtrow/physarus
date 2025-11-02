use crate::prelude::*;

// Step 1: Multiple Importance Sampling Weight Derivation
// In MIS, we combine estimates from multiple sampling strategies to reduce variance.
// For two strategies f and g with nf and ng samples respectively,
// we need weights wf and wg such that wf + wg = 1.
//
// Step 2: Balance Heuristic Formulation
// The balance heuristic assigns weights proportional to the "strength" of each strategy,
// where strength = (number of samples) × (probability density function value).
//
// For strategy f: strength_f = nf × pdf_f
// For strategy g: strength_g = ng × pdf_g
//
// Step 3: Normalization to Ensure wf + wg = 1
// wf = strength_f / (strength_f + strength_g)
// wf = (nf × pdf_f) / (nf × pdf_f + ng × pdf_g)
//
// This is the balance heuristic weight for strategy f.

/// Computes balance heuristic weight for strategy f in multiple importance sampling.
/// Used to optimally combine samples from two different sampling strategies.
///
/// Parameters: `nf` (samples from f), `pdf_f` (density of f), `ng` (samples from g), `pdf_g` (density of g)
/// Returns weight for strategy `f: nf*pdf_f / (nf*pdf_f + ng*pdf_g)`
pub fn balance_heuristic(nf: usize, pdf_f: Probpart, ng: usize, pdf_g: Probpart) -> Float {
    let (nf, ng) = (nf as Float, ng as Float);

    nf * pdf_f / (nf * pdf_f + ng * pdf_g)
}

// Step 1: Power Heuristic with β = 2
// The power heuristic is a generalization where weights use powers of the sampling densities.
// For power β and strategy f: wf = (nf × pdf_f)^β / Σ(ni × pdf_i)^β
//
// The power heuristic with β = 2 often provides better variance reduction
// than the balance heuristic when one strategy is significantly better.

/// Computes power heuristic weight (β=2) for strategy f in multiple importance sampling.  
/// Generally provides better variance reduction than balance heuristic when strategies differ significantly.
///
/// Parameters: `nf` (samples from f), `pdf_f` (density of f), `ng` (samples from g), `pdf_g` (density of g)  
/// Returns weight for strategy `f: (nf*pdf_f)² / ((nf*pdf_f)² + (ng*pdf_g)²)`
pub fn power_heuristic(nf: usize, pdf_f: Probpart, ng: usize, pdf_g: Probpart) -> Float {
    let (nf, ng) = (nf as Float, ng as Float);

    let f = nf * pdf_f;
    let g = ng * pdf_g;

    f * f / (f * f + g * g)
}

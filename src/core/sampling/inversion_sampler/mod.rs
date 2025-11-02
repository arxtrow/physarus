//! # Inversion Sampling Algorithms
//!
//! This module implements multiple inversion-sampling techniques for
//! discrete and continuous probability distributions. See submodules
//! for documentation and theory proofs for each specific distribution.
//!
//! General idea for Discrete and Continuous case
//! ```text
//! Discrete Inversion (4 Steps)
//! 1. Start with unnormalized weights w[i] ≥ 0 for indices i = 0, 1, ..., n-1
//! 3. Build discrete CDF: CDF[i] = Σ(j=0 to i) w[j] = cumulative probability up to index i
//! 4. Sample: u ~ Uniform(0,1), find i such that CDF[i-1] < u ≤ CDF[i], return i
//!
//! Why it works: If U ~ Uniform(0,1) and i = CDF⁻¹(U)
//! then i is distributed according to PMF p[i] = w[i] / weights_sum
//!
//! Proof: Need to show P(I = i) = p[i] for all indices i
//! By construction: P(I = i) = P(CDF[i-1] < U ≤ CDF[i])    // Definition of discrete inversion
//!                           = CDF[i] - CDF[i-1]           // U ~ Uniform(0,1): P(a < U ≤ b) = b - a
//!                           = p[i]                        // Definition of discrete CDF
//!
//! Since P(I = i) = p[i] for all i, the sampled index I follows the discrete distribution with PMF p[i]
//! ```
//! ```text
//! Continuous Inversion (4 Steps)
//! 1. Start with unnormalized density f(x) ≥ 0 over domain D
//! 2. Normalize: p(x) = f(x)/c where c = ∫D f(x) dx
//! 3. Build CDF: CDF(x) = ∫{x₀}^x p(t) dt = (1/c) ∫{x₀}^x f(t) dt
//! 4. Sample: u ~ Uniform(0,1), solve x = CDF⁻¹(u)
//!
//!
//! Why it works
//!
//! If U ~ Uniform(0,1) and X = CDF⁻¹(U), then X follows the distribution with PDF p(X)
//!
//! Proof:
//! Given: U ~ Uniform(0,1) and X = CDF⁻¹(U)
//! What We Want to Prove: X = CDF⁻¹(U) has the same distribution as a random variable with PDF p(x), Y ~ p(x)
//!
//! Prerequisites:
//! - f(x) is integrable over D (so normalization constant c exists and is finite)
//! - CDF is strictly increasing (so CDF⁻¹ exists and is unique)
//!
//! Proof:
//! By definition, if Y ~ p(x), then: P(Y ≤ x) = ∫_{-∞}^x p(t) dt = CDF(x)
//!
//! P(X ≤ x) = P(CDF⁻¹(U) ≤ x)    // Definition of X
//!          = P(U ≤ CDF(x))      // Apply CDF to both sides (CDF monotonic: x₁ ≤ x₂ => f(x₁) ≤ f(x₂))
//!          = CDF(x)             // Since U ~ Uniform(0,1), P(U ≤ y) = y for y ∈ [0,1]
//!
//! Therefore: P(X ≤ x) = CDF(x) = P(Y ≤ x), X and Y have identical CDFs => they have identical distributions => X ~ p(x)
//! ```
pub mod bilinear;
pub mod discrete;
pub mod disk;
pub mod exponential;
pub mod linear;
pub mod logistic;
pub mod normal;
pub mod smooth_step;
pub mod tent;

#[derive(Debug, Clone)]
pub struct InversionSampler;

// These use their own structs instead of InversionSampler
pub mod piecewise_const_1d;

use std::collections::VecDeque;

use crate::prelude::*;

/// AliasSampler: O(1) discrete sampling from arbitrary weights distributions.
///
/// Properties
/// - Preprocessing is O(n)
/// - Each sample is O(1), unbiased, and matches the original distribution.
/// - Remapping ensures correct downstream stratification for compositional/inverse sampling.
#[derive(Default, Debug, Clone)]
struct AliasSampler {
    bins: Vec<Bin>,
}

#[derive(Default, Debug, Clone)]
/// One probability bucket for alias sampling.
/// - `prob`: PMF of the original outcome
/// - `threshold`: split point between `prob` and its alias's
/// - `alias_idx`: other bucket used when threshold surpassed
struct Bin {
    prob: Probpart,
    threshold: Probpart,
    alias_idx: Option<usize>,
}

impl AliasSampler {
    pub fn new(weights: &[Float]) -> Self {
        let weights_sum: Float = weights.iter().sum();

        let mut bins: Vec<Bin> = weights
            .iter()
            .map(|&weight| Bin {
                prob: Probpart::new((weight as f64 / weights_sum as f64) as Float), // normalize weights to probabilities
                threshold: Probpart::ONE,
                alias_idx: None,
            })
            .collect();

        let bins_len = bins.len();

        // Rescale prob so avg(scaled_probs) == 1.0
        let scaled_probs: Vec<Float> = bins
            .iter()
            .map(|bin| bin.prob * bins_len as Float)
            .collect();

        // Partitioned by (<1.0) and (>=1.0)
        // note that indeces taken from `scaled_props`
        let (mut underfull, mut overfull): (VecDeque<_>, VecDeque<_>) =
            scaled_probs.iter().enumerate().fold(
                (VecDeque::new(), VecDeque::new()),
                |(mut left, mut right), (i, &sp)| {
                    if sp < 1.0 {
                        left.push_back((i, sp));
                    } else {
                        right.push_back((i, sp))
                    }

                    (left, right)
                },
            );

        while let (Some((u_idx, u_p)), Some((o_idx, o_p))) =
            (underfull.pop_front(), overfull.pop_front())
        {
            // Assign threshold and alias so each bin contains its own probability up to 1.
            bins[u_idx].threshold = Probpart::new(u_p);
            bins[u_idx].alias_idx = Some(o_idx);

            let p_excess = o_p + u_p - 1.0;
            if p_excess >= 1.0 {
                overfull.push_back((o_idx, p_excess))
            } else {
                underfull.push_back((o_idx, p_excess))
            }
        }

        // Any remaining bins exist due to floating-point round-off error;
        // These items have `p` slightly less or slightly greater than 1.0
        // so they should be given threshold = 1.0 with no alias
        for (i, _) in underfull.into_iter().chain(overfull.into_iter()) {
            bins[i].threshold = Probpart::ONE;
            bins[i].alias_idx = None;
        }

        Self { bins }
    }

    /// Returns (`idx`, `pmf`, `u_remapped`)
    /// `idx`: index of sample from original weights distribution
    /// `pmf`: probability of chosen index: weights[idx] / weights_sum
    /// `u_remapped`: a new uniform random variable in [0, 1)
    ///
    pub fn sample(&self, u: Probpart) -> (usize, Probpart, Probpart) {
        let n = self.bins.len() as Float;

        let bin_idx = (u * n).min(n - 1.0).floor() as usize;

        let u_part_in_bin = Probpart::new(((u * n) - bin_idx as Float).min(CLOSEST_TO_ONE));

        let bin = &self.bins[bin_idx];

        if u_part_in_bin < bin.threshold {
            (bin_idx, bin.prob, u_part_in_bin / bin.threshold)
        } else {
            let alias_idx = bin.alias_idx.expect(
                "by construction in AliasSampler::new, alias_idx is None only if bin.threshold = 1.0
                 and we cannot be in this branch if bin.threshold = 1.0",
            );

            let alias_bin = &self.bins[alias_idx];

            (
                alias_idx,
                alias_bin.prob,
                (u_part_in_bin - bin.threshold) / (Probpart::ONE - bin.threshold),
            )
        }
    }

    pub fn pmf(&self, idx: usize) -> Option<Probpart> {
        self.bins.get(idx).and_then(|bin| Some(bin.prob))
    }
}

// #[cfg(feature = "samplers_tests")]
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn empirical_distribution_matches_input_probabilities() {
        let weights = vec![5.5, 9.2, 10.02, 0.1, 2.88, 1.05, 1.05];
        let sampler = AliasSampler::new(&weights);

        const N_SAMPLES: usize = 1_000_000;
        let mut counts = vec![0; weights.len()];

        let mut rng = rand::rng();

        for _ in 0..N_SAMPLES {
            let u: Float = rng.random();
            let (idx, _, _) = sampler.sample(Probpart::new(u));
            counts[idx] += 1;
        }

        let total: Float = weights.iter().sum();
        let expected_probs: Vec<Float> = weights.iter().map(|&w| w / total).collect();

        let empirical: Vec<Float> = counts
            .iter()
            .map(|&c| c as Float / N_SAMPLES as Float)
            .collect();

        let tolerance = 0.002;

        for (i, (exp, obs)) in expected_probs.iter().zip(&empirical).enumerate() {
            assert!(
                (obs - exp).abs() < tolerance,
                "Prob mismatch for idx {}: expected {:.4}, got {:.4}",
                i,
                exp,
                obs
            );
        }
    }
}

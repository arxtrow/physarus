use crate::core::{CLOSEST_TO_ONE, Float, Probpart};

impl super::InversionSampler {
    /// Given uniform random variable `u` in `[0, 1)`
    /// samples an index from discrete distribution defined by `weights`
    /// Over many calls with independent uniform `u`,
    /// returned indices are distributed proportionally to the weights.
    ///
    /// Returns (`idx`, `pmf`, `u_remapped`)
    /// `idx`: index sampled by mapping `u` onto the weights distribution
    /// `pmf`: probability of chosen index: `weights[idx] / weights_sum`
    /// `u_remapped`: a new uniform random variable in `[0, 1)`
    ///
    pub fn sample_discrete(weights: &[Float], u: Probpart) -> (usize, Probpart, Probpart) {
        let weights_sum: Float = weights.iter().sum();
        assert!(weights_sum >= 0.0);

        let target = (u * weights_sum).min(weights_sum.next_down());

        let cdf = weights.iter().scan(0.0, |sum, &w| {
            *sum += w;
            Some(*sum)
        });

        let (idx, sum) = cdf
            .enumerate()
            .find(|(_, sum)| *sum > target)
            .expect("target <= weights_sum.next_down() < (last prefix_sum)");

        let weight = weights[idx];
        let pmf = weight / weights_sum;
        let prev_sum = sum - weight;
        let u_remapped = ((target - prev_sum) / weight).min(CLOSEST_TO_ONE);

        return (idx, pmf.into(), u_remapped.into());
    }
}

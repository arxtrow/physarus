use crate::prelude::*;

#[derive(Debug, Clone)]
pub struct ReservoirSampler<T: Clone> {
    weight_sum: Float,
    reservoir_weight: Float,
    reservoir: Option<T>,
}

impl<T: Clone> ReservoirSampler<T> {
    pub fn new() -> Self {
        Self {
            weight_sum: 0.0,
            reservoir_weight: 0.0,
            reservoir: None,
        }
    }

    pub fn add(&mut self, sample: T, weight: Float) {
        self.weight_sum += weight;

        let p = weight / self.weight_sum;

        if rand::random::<Float>() < p {
            self.reservoir = Some(sample);
            self.reservoir_weight = weight;
        }
    }

    pub fn add_with(&mut self, cb: fn() -> T, weight: Float) {
        self.weight_sum += weight;

        let p = weight / self.weight_sum;

        if rand::random::<Float>() < p {
            self.reservoir = Some(cb());
            self.reservoir_weight = weight;
        }
    }

    pub fn has_sample(&self) -> bool {
        self.weight_sum > 0.0 && self.reservoir.is_some()
    }

    pub fn get_sample(&self) -> &Option<T> {
        &self.reservoir
    }

    pub fn sample_prob(&self) -> Probpart {
        Probpart::new(self.reservoir_weight / self.weight_sum)
    }

    pub fn reset(&mut self) {
        self.weight_sum = 0.0;
        self.reservoir_weight = 0.0;
    }

    pub fn merge(&mut self, rs: &Self) {
        if rs.has_sample() {
            self.add(rs.reservoir.clone().expect("checked in if"), rs.weight_sum);
        }
    }
}

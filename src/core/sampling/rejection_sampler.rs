use crate::prelude::*;

#[derive(Debug, Clone)]
pub struct RejectionSampler;

impl RejectionSampler {
    pub fn sample_disk() -> Point2 {
        let (x, y) = loop {
            let x = 1.0 - 2.0 * rand::random::<Float>();
            let y = 1.0 - 2.0 * rand::random::<Float>();

            if x * x + y * y <= 1.0 {
                break (x, y);
            }
        };

        Point2::new(x, y)
    }
}

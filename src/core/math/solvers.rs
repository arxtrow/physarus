use crate::prelude::*;

// Step 1: Standard Quadratic Equation
// ax² + bx + c = 0, where a ≠ 0
//
// Step 2: Complete the Square Method
// ax² + bx + c = 0
// a(x² + (b/a)x) + c = 0
// a(x² + (b/a)x + (b/2a)²) - a(b/2a)² + c = 0
// a(x + b/2a)² = b²/4a - c
// (x + b/2a)² = (b² - 4ac)/(4a²)
//
// Step 3: Derive Standard Quadratic Formula
// x + b/2a = ±√((b² - 4ac)/(4a²))
// x = -b/2a ± √(b² - 4ac)/(2a)
// x = (-b ± √(b² - 4ac))/(2a)
//
// Step 4: Discriminant Analysis
// Δ = b² - 4ac
// If Δ < 0: no real solutions
// If Δ = 0: one repeated real solution
// If Δ > 0: two distinct real solutions
//
// Step 5: Numerical Stability Problem
// When b and √Δ have same sign, computing -b ± √Δ causes catastrophic cancellation
// Example: If b = 1000, √Δ = 999.999, -b + √Δ = -0.001, but both inputs
// might carry tiny rounding errors, yielding a large relative error in the result.
//
// Step 6: Citardauq Formula (Numerically Stable Alternative)
// Define: q = -½(b + sign(b)√Δ)
// This ensures we always subtract terms of same sign, avoiding cancellation
//
// From Vieta's formulas: x₁x₂ = c/a and x₁ + x₂ = -b/a
// If x₁ = q/a, then x₂ = (c/a)/x₁ = c/q
// Final formulas: x₁ = q/a, x₂ = c/q

/// Solves quadratic equation ax² + bx + c = 0 using numerically stable methods.
/// Handles degenerate cases when a=0 (linear equation) or a=b=0 (no solution).
///
/// Returns `Some((t0, t1))` where t0 ≤ t1 are the real roots, or `None` if no real solutions exist.
/// Uses Citardauq formula to avoid catastrophic cancellation in floating-point arithmetic.
#[inline]
pub fn quadratic(a: Float, b: Float, c: Float) -> Option<(Float, Float)> {
    let (a, b, c) = (a as f64, b as f64, c as f64);

    match (a, b, c) {
        (0.0, 0.0, _) => None,

        (0.0, b, c) => Some(((-c / b) as Float, (-c / b) as Float)),

        (a, b, c) => {
            let discrim = (b * b - 4.0 * a * c).ensure(|&d| d >= 0.0, None)?;
            let root_discrim = discrim.sqrt();

            let q = -0.5 * (b + b.signum() * root_discrim);

            let (t0, t1) = match (q / a, c / q) {
                (l, r) if l > r => (r, l),
                (l, r) => (l, r),
            };

            Some((t0 as Float, t1 as Float))
        }
    }
}

// Step 1: Root Finding Problem Setup
// We seek to solve f(x) = 0 for a continuous function f on interval [x_left, x_right]
// where f(x_left) and f(x_right) have opposite signs.
// By the Intermediate Value Theorem, a root exists in the interval.
//
// Step 2: Secant Method Initial Estimate
// The secant line through points (x_left, f_left) and (x_right, f_right) has slope:
// m = (f_right - f_left) / (x_right - x_left)
// and secant line equation
// s(x) = f_left + m(x - x_left)
//
// Setting the secant line equal to zero to find x-intercept:
// 0 = f_left + m(x - x_left)
// x = x_left - f_left/m
// x = x_left - f_left(x_right - x_left)/(f_right - f_left)
// or after rearranging
// x = (x_left·f_right - f_left·x_right) / (f_right - f_left) this value is barycentric balance of the f_left and f_right
//
// Step 3: Newton's Method Iteration
// Newton's method uses linear approximation around current point x₀:
// f(x) ≈ f(x₀) + f'(x₀)(x - x₀)
// Setting f(x) = 0 and solving for x:
// 0 = f(x₀) + f'(x₀)(x - x₀)
// x = x₀ - f(x₀)/f'(x₀)
//
// Step 4: Bisection Method Bracket Maintenance
// To ensure convergence, maintain bracket [x_left, x_right] where f(x_left).signum() != f(x_right).signum()
// Update rule based on sign of f(x_mid):
// If sign(f(x_mid)) = sign(f(x_left)): x_left = x_mid
// If sign(f(x_mid)) = sign(f(x_right)): x_right = x_mid
//
// Key formulas in implementation:
// Initial estimate: x_mid = (f_right·x_left - f_left·x_right)/(f_right - f_left)
// Bisection fallback: x_mid = (x_left + x_right)/2
// Newton update: x_mid = x_mid - f(x_mid)/f'(x_mid)

/// Finds root of function using hybrid Newton-bisection method.
/// Combines Newton's quadratic convergence with bisection's guaranteed convergence.
///
/// `f`: function returning (value, derivative) at given point
/// `x_left`, `x_right`: initial bracket where f(x_left) and f(x_right) have opposite signs
///
/// Returns root x where f(x) ≈ 0
#[inline]
pub fn newton_bisection<F>(f: F, mut x_left: Float, mut x_right: Float) -> Float
where
    F: Fn(Float) -> (Float, Float), // Returns (function_value, derivative)
{
    assert!(x_left < x_right);

    // Get and Check function endpoints for roots
    let (f_left, f_right) = {
        let (f_left, _) = f(x_left);
        if f_left.is_zero_eps() {
            return x_left;
        }

        let (f_right, _) = f(x_right);
        if f_right.is_zero_eps() {
            return x_right;
        }

        assert!(
            f_left.signum() != f_right.signum(),
            "Function values at bounds must have opposite signs for a guaranteed root"
        );

        (f_left, f_right)
    };

    // Initial guess using secant method for faster start
    let mut x_mid = (f_right * x_left - f_left * x_right) / (f_right - f_left);

    loop {
        if !x_mid.is_finite() || !(x_left..x_right).contains(&x_mid) {
            x_mid = (x_left + x_right) / 2.0;
        }

        let (fx_mid, dfx_mid) = f(x_mid);

        // Bisection: update interval based on function signs
        if f_left.signum() == fx_mid.signum() {
            x_left = x_mid;
        } else {
            x_right = x_mid
        }

        // Check convergence
        if x_right.eq_eps(&x_left) || fx_mid.abs().is_zero_eps() {
            return x_mid;
        }

        // Newton step
        x_mid = x_mid - fx_mid / dfx_mid;
    }
}

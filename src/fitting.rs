//! Logistic regression fitting for duplication rate vs expression.
//!
//! Implements the GLM binomial logit fit that dupRadar uses:
//!   logit(dupRate) = β₀ + β₁ * log10(RPK)
//!
//! Uses Iteratively Reweighted Least Squares (IRLS) to fit the model,
//! matching R's `glm(y ~ x, family = binomial(link = "logit"))`.

use anyhow::Result;

/// Result of the logistic regression fit.
#[derive(Debug, Clone)]
pub struct FitResult {
    /// Raw intercept coefficient (β₀)
    pub beta0: f64,
    /// Raw slope coefficient (β₁)
    pub beta1: f64,
    /// exp(β₀) - the value dupRadar reports as "intercept"
    pub intercept: f64,
    /// exp(β₁) - the value dupRadar reports as "slope"
    pub slope: f64,
}

impl FitResult {
    /// Predict duplication rate for a given log10(RPK) value.
    ///
    /// dupRate = 1 / (1 + exp(-(β₀ + β₁ * x)))
    pub fn predict(&self, log10_rpk: f64) -> f64 {
        let eta = self.beta0 + self.beta1 * log10_rpk;
        1.0 / (1.0 + (-eta).exp())
    }

    /// Predict duplication rate for raw RPK value.
    pub fn predict_rpk(&self, rpk: f64) -> f64 {
        if rpk <= 0.0 {
            return self.predict(f64::NEG_INFINITY);
        }
        self.predict(rpk.log10())
    }
}

/// Fit a logistic regression model: dupRate ~ log10(RPK)
///
/// Uses IRLS (Iteratively Reweighted Least Squares) to solve the GLM
/// with binomial family and logit link, matching R's `glm()` behavior.
///
/// # Arguments
/// * `rpk` - Per-gene RPK values
/// * `dup_rate` - Per-gene duplication rates (0.0 to 1.0, may contain NaN)
///
/// # Returns
/// FitResult with coefficients matching dupRadar's `duprateExpFit()` output.
pub fn duprate_exp_fit(rpk: &[f64], dup_rate: &[f64]) -> Result<FitResult> {
    anyhow::ensure!(
        rpk.len() == dup_rate.len(),
        "rpk and dup_rate must have the same length ({} vs {})",
        rpk.len(),
        dup_rate.len()
    );
    // Filter to genes with valid data: RPK > 0 and finite dupRate
    let mut x_vals: Vec<f64> = Vec::new();
    let mut y_vals: Vec<f64> = Vec::new();

    for (r, d) in rpk.iter().zip(dup_rate.iter()) {
        if *r > 0.0 && d.is_finite() {
            x_vals.push(r.log10());
            // Clamp y to avoid exact 0 and 1 (which cause issues with logit)
            let y_clamped = d.clamp(1e-10, 1.0 - 1e-10);
            y_vals.push(y_clamped);
        }
    }

    let n = x_vals.len();
    anyhow::ensure!(n >= 2, "Need at least 2 valid data points for fitting");

    // IRLS for logistic regression
    // Model: logit(μ) = β₀ + β₁ * x
    // where logit(μ) = ln(μ/(1-μ))
    //
    // IRLS iteration:
    //   1. Compute linear predictor: η = Xβ
    //   2. Compute mean: μ = logistic(η)
    //   3. Compute variance: v = μ(1-μ)
    //   4. Compute working response: z = η + (y-μ)/v
    //   5. Compute weights: w = v
    //   6. Solve weighted least squares: β_new = (X'WX)^{-1} X'Wz

    // Match R's glm() defaults for convergence
    let max_iter: usize = 25; // R default: control$maxit = 25
    let tol: f64 = 1e-8; // R default: control$epsilon = 1e-8

    // Initialize with starting values (R's default initialization)
    // Start with the mean of y transformed through the link
    let y_mean: f64 = y_vals.iter().sum::<f64>() / n as f64;
    let y_mean = y_mean.clamp(0.01, 0.99);
    let eta_start = (y_mean / (1.0 - y_mean)).ln();

    let mut beta0 = eta_start;
    let mut beta1 = 0.0;

    let mut converged = false;
    for _iter in 0..max_iter {
        // Compute η, μ, and working quantities
        let mut sw = 0.0; // sum of weights
        let mut swx = 0.0; // sum of w*x
        let mut swx2 = 0.0; // sum of w*x^2
        let mut swz = 0.0; // sum of w*z
        let mut swxz = 0.0; // sum of w*x*z

        for i in 0..n {
            let xi = x_vals[i];
            let yi = y_vals[i];

            // Linear predictor
            let eta = beta0 + beta1 * xi;

            // Mean (inverse logit)
            let mu = 1.0 / (1.0 + (-eta).exp());

            // Variance function for binomial
            let v = mu * (1.0 - mu);

            // Guard against near-zero variance
            if v < 1e-20 {
                continue;
            }

            // Weight
            let w = v;

            // Working response
            let z = eta + (yi - mu) / v;

            sw += w;
            swx += w * xi;
            swx2 += w * xi * xi;
            swz += w * z;
            swxz += w * xi * z;
        }

        // Solve 2x2 weighted normal equations:
        // [sw    swx ] [β₀]   [swz ]
        // [swx   swx2] [β₁] = [swxz]
        let det = sw * swx2 - swx * swx;
        if det.abs() < 1e-30 {
            anyhow::bail!("Singular matrix in IRLS iteration");
        }

        let new_beta0 = (swx2 * swz - swx * swxz) / det;
        let new_beta1 = (sw * swxz - swx * swz) / det;

        // Check convergence
        let d0 = (new_beta0 - beta0).abs();
        let d1 = (new_beta1 - beta1).abs();

        beta0 = new_beta0;
        beta1 = new_beta1;

        if d0 < tol && d1 < tol {
            converged = true;
            break;
        }
    }

    if !converged {
        log::warn!(
            "IRLS logistic regression did not converge within {} iterations (tolerance: {:.0e})",
            max_iter,
            tol
        );
    }

    Ok(FitResult {
        beta0,
        beta1,
        intercept: beta0.exp(),
        slope: beta1.exp(),
    })
}

/// Compute the RPKM threshold line position in log10(RPK) space.
///
/// Matches dupRadar's calculation: find the RPK value at which RPKM ≈ threshold.
/// Due to a quirk in dupRadar's R code, this effectively returns the lowest RPK
/// value where RPKM >= threshold.
pub fn compute_rpkm_threshold_rpk(rpk: &[f64], rpkm: &[f64], threshold: f64) -> Option<f64> {
    // Find genes with RPKM >= threshold, get the one with lowest RPK
    let mut rpk_gt: Option<f64> = None;
    for (r, m) in rpk.iter().zip(rpkm.iter()) {
        if *m >= threshold && *r > 0.0 {
            match rpk_gt {
                None => rpk_gt = Some(*r),
                Some(current) => {
                    if *r < current {
                        rpk_gt = Some(*r);
                    }
                }
            }
        }
    }

    rpk_gt
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_predict() {
        let fit = FitResult {
            beta0: -2.0,
            beta1: 1.5,
            intercept: (-2.0f64).exp(),
            slope: 1.5f64.exp(),
        };

        // At x=0 (RPK=1): p = 1/(1+exp(2)) ≈ 0.119
        let p = fit.predict(0.0);
        assert!((p - 0.1192).abs() < 0.01);

        // At very high x, should approach 1.0
        let p = fit.predict(100.0);
        assert!(p > 0.99);
    }

    #[test]
    fn test_fit_basic() {
        // Create synthetic data: clear logistic relationship
        let n = 1000;
        let mut rpk = Vec::with_capacity(n);
        let mut dup_rate = Vec::with_capacity(n);

        // True model: logit(p) = -3 + 1.5 * log10(rpk)
        for i in 0..n {
            let log_rpk = -1.0 + 5.0 * (i as f64 / n as f64);
            let r = 10.0f64.powf(log_rpk);
            let eta = -3.0 + 1.5 * log_rpk;
            let p = 1.0 / (1.0 + (-eta).exp());

            rpk.push(r);
            dup_rate.push(p);
        }

        let fit = duprate_exp_fit(&rpk, &dup_rate).unwrap();

        // Should recover approximately β₀=-3, β₁=1.5
        assert!(
            (fit.beta0 - (-3.0)).abs() < 0.1,
            "beta0={}, expected -3.0",
            fit.beta0
        );
        assert!(
            (fit.beta1 - 1.5).abs() < 0.1,
            "beta1={}, expected 1.5",
            fit.beta1
        );
    }

    #[test]
    fn test_fit_with_nans() {
        // Should handle NaN dupRate values (genes with 0 reads)
        let rpk = vec![10.0, 100.0, 1000.0, 0.0, 50.0];
        let dup_rate = vec![0.1, 0.3, 0.7, f64::NAN, 0.2];

        let fit = duprate_exp_fit(&rpk, &dup_rate);
        assert!(fit.is_ok());
    }

    #[test]
    fn test_rpkm_threshold() {
        let rpk = vec![10.0, 50.0, 100.0, 500.0, 1000.0];
        let rpkm = vec![0.1, 0.4, 0.6, 2.0, 5.0];

        let thresh = compute_rpkm_threshold_rpk(&rpk, &rpkm, 0.5);
        // The lowest RPK where RPKM >= 0.5 is RPK=100 (RPKM=0.6)
        assert_eq!(thresh, Some(100.0));
    }
}

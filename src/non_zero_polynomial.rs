//! Module for binary polynomial that is guaranteed to be non-zero

use std::fmt::Display;
use num_traits::{One, Zero};
use crate::BinaryPolynomial;

/// Represents a non-zero binary univariate polynomial
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct NonZeroBinaryPolynomial(BinaryPolynomial);

impl NonZeroBinaryPolynomial {

    /// Creates a new `NonZeroBinaryPolynomial` from a binary polynomial
    ///
    /// # Arguments
    /// * `value` - The binary polynomial value to create a `NonZeroBinaryPolynomial` from
    ///
    /// # Returns
    /// `None` if the value is zero, `Some(NonZeroBinaryPolynomial)` containing the value of the passed polynomial otherwise
    pub fn new(value: BinaryPolynomial) -> Option<Self> {
        if value.is_zero() {
            None
        } else {
            Some(NonZeroBinaryPolynomial(value))
        }
    }

    /// Returns the binary polynomial value of the non-zero polynomial
    ///
    /// # Returns
    /// The binary polynomial value of the non-zero polynomial
    pub fn get(&self) -> &BinaryPolynomial {
        &self.0
    }

    /// Returns the owned binary polynomial value of the non-zero polynomial
    ///
    /// # Returns
    /// The owned binary polynomial value of the non-zero polynomial
    pub fn get_owned(self) -> BinaryPolynomial {
        self.0
    }

    /// Check if the non-zero polynomial is coprime with the other non-zero polynomial
    ///
    /// # Arguments
    /// * `other` - The other non-zero polynomial to check coprimality with
    ///
    /// # Returns
    /// `true` if the non-zero polynomials are coprime, `false` otherwise
    pub fn coprime(&self, other: &Self) -> bool {
        self.gcd(other).get().is_one()
    }

    /// Returns the Greatest Common Divisor of the non-zero polynomials
    ///
    /// # Arguments
    /// * `other` - The other non-zero polynomial to calculate the GCD with
    ///
    /// # Returns
    /// The Greatest Common Divisor of the non-zero polynomials (which is always non-zero)
    pub fn gcd(&self, other: &Self) -> Self {
        let mut a = self.get().clone();
        let mut b = other.get().clone();

        while let Some(non_zero_b) = Self::new(b.clone()) {
            (a, b) = (b, a % non_zero_b);
        }

        Self::new(a).unwrap()
    }

    /// Computes the Extended Greatest Common Divisor of two non-zero polynomials
    ///
    /// # Arguments
    /// * `other` - The other non-zero polynomial to calculate the Extended Greatest Common Divisor with
    ///
    /// # Returns
    /// A tuple `(d, x, y)` where `d` is the Greatest Common Divisor of polynomials `self` and `other`, `x`, `y` are polynomials that satisfy: `(self * x) + (other * y) = d`
    pub fn egcd(&self, other: &Self) -> (Self, BinaryPolynomial, BinaryPolynomial) {
        let mut a: (Self, BinaryPolynomial, BinaryPolynomial) = (self.clone(), BinaryPolynomial::one(), BinaryPolynomial::zero());
        let mut b: (Self, BinaryPolynomial, BinaryPolynomial) = (other.clone(), BinaryPolynomial::zero(), BinaryPolynomial::one());

        loop {
            let (q, r) = a.0.get().div_mod(&b.0);
            if r.is_zero() {
                return b;
            }
            (a, b) = (b.clone(),
                      (NonZeroBinaryPolynomial::new(r).unwrap(), BinaryPolynomial::from(a.1.polynomial ^ (q.clone() * b.1).polynomial), BinaryPolynomial::from(a.2.polynomial ^ (q * b.2).polynomial)));
        }
    }

    /// Computes the modular inverse of the non-zero polynomial
    ///
    /// # Arguments
    /// * `modulus` - The non-zero modulus polynomial
    ///
    /// # Returns
    /// The modular inverse of the non-zero polynomial, or `None` if the polynomial is not invertible modulo the modulus polynomial
    pub fn inv_mod(&self, modulus: &Self) -> Option<BinaryPolynomial> {
        let (d, x, _) = self.egcd(modulus);
        if !d.get().is_one() {
            return None;
        }
        Some(x)
    }
}

/// Display implementation for `NonZeroBinaryPolynomial`
///
/// The polynomial is displayed in the form of a sum of monomials, like: `x^3 + x^2 + 1`, or `0` for the zero polynomial
///
/// This is the same as the display implementation for `BinaryPolynomial`
impl Display for NonZeroBinaryPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};
    use crate::{BinaryPolynomial, NonZeroBinaryPolynomial};

    #[test]
    fn test_gcd() {
        let polynomial = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, true, true, false, true])).unwrap();
        let polynomial2 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, true])).unwrap();

        let result = polynomial.gcd(&polynomial2);
        assert_eq!(result.to_string(), "x + 1");

        let result = polynomial2.gcd(&polynomial);
        assert_eq!(result.to_string(), "x + 1");

        let result = polynomial.gcd(&polynomial);
        assert_eq!(result, polynomial);
    }

    #[test]
    fn test_egcd() {
        let polynomial = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, true, true, false, true])).unwrap();
        let polynomial2 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, true])).unwrap();

        let result = polynomial.egcd(&polynomial2);
        assert_eq!(result.0.to_string(), "x + 1");
        assert_eq!(result.1.to_string(), "1");
        assert_eq!(result.2.to_string(), "x^2 + x");

        let result = polynomial2.egcd(&polynomial);
        assert_eq!(result.0.to_string(), "x + 1");
        assert_eq!(result.1.to_string(), "x^2 + x");
        assert_eq!(result.2.to_string(), "1");

        let result = polynomial.egcd(&polynomial);
        assert_eq!(result.0.to_string(), "x^4 + x^3 + x^2 + 1");
        assert_eq!(result.1.to_string(), "0");
        assert_eq!(result.2.to_string(), "1");
    }

    #[test]
    fn test_inv_mod() {
        let modulus_polynomial = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, false, true])).unwrap();
        let polynomial = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, true, true])).unwrap();
        let polynomial2 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, true])).unwrap();
        let polynomial3 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, false])).unwrap();

        let result = polynomial.inv_mod(&modulus_polynomial).unwrap();
        assert_eq!(result.to_string(), "x^3 + x + 1");

        let result = polynomial2.inv_mod(&modulus_polynomial);
        assert!(result.is_none());

        let result = polynomial3.inv_mod(&modulus_polynomial).unwrap();
        assert_eq!(result.to_string(), "x");
    }

    #[test]
    fn test_nonzeropolynomial() {
        let zero_polynomial = BinaryPolynomial::from(vec![false, false, false, false]);
        assert_eq!(NonZeroBinaryPolynomial::new(zero_polynomial), None);

        let polynomial = BinaryPolynomial::zero();
        assert_eq!(NonZeroBinaryPolynomial::new(polynomial), None);

        let polynomial = BinaryPolynomial::one();
        assert_eq!(NonZeroBinaryPolynomial::new(polynomial.clone()).unwrap().get(), &polynomial);

        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
        assert_eq!(NonZeroBinaryPolynomial::new(polynomial.clone()).unwrap().get(), &polynomial);
    }

    #[test]
    fn test_coprime() {
        let polynomial = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![false, true, false, true])).unwrap();
        let polynomial2 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, true])).unwrap();
        let polynomial3 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, true, false, true])).unwrap();

        let result = polynomial.coprime(&polynomial2);
        assert!(!result);

        let result = polynomial.coprime(&polynomial);
        assert!(!result);

        let result = polynomial2.coprime(&polynomial3);
        assert!(result);
    }
}
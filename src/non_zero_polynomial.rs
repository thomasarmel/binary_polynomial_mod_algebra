//! Module for binary polynomial that is guaranteed to be non-zero

use crate::berkelamp_matrix::BerkelampMatrix;
use crate::error::BinaryPolynomialError;
use crate::BinaryPolynomial;
use num_traits::{One, Zero};
use std::collections::{HashMap, HashSet, VecDeque};
use std::fmt::Display;
use std::ops::{AddAssign, Deref, Mul};
use num_bigint::BigUint;
use crate::utils::prime_factors;

/// Represents a non-zero binary univariate polynomial
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    /// Returns the binary polynomial value of the non-zero polynomial
    ///
    /// # Returns
    /// The binary polynomial value of the non-zero polynomial
    pub fn get_mut(&mut self) -> &mut BinaryPolynomial {
        &mut self.0
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

    /// Check if the polynomial is irreducible using Rabin's irreducibility test
    pub fn is_irreducible(&self) -> bool {
        let degree = self.degree();
        if self.is_one() {
            return false;
        }
        let degree = degree as usize;
        let x_polynomial = BinaryPolynomial::from(BigUint::from(2usize)); // P(X) = X
        let self_nonzero = self.clone();

        for q in prime_factors(degree) {
            let h = x_polynomial.pow_mod(1 << (degree / q), &self_nonzero) + (x_polynomial.clone() % self_nonzero.clone());
            let h_nonzero = NonZeroBinaryPolynomial::new(h).unwrap();
            if self_nonzero.gcd(&h_nonzero) != Self::one() {
                return false
            }
        }
        let h = x_polynomial.pow_mod(1 << degree, &self_nonzero) + (x_polynomial % self_nonzero);
        h.is_zero()
    }

    /// Check if the polynomial is primitive
    ///
    /// A binary polynomial of degree m is primitive if it is irreducible
    /// and if the smallest positive integer n such that the polynomial divides X^n + 1 is n = 2^m - 1.
    pub fn is_primitive(&self) -> bool {
        if !self.is_irreducible() {
            return false;
        }
        let x_polynomial = BinaryPolynomial::from(BigUint::from(2usize)); // P(X) = X
        if **self == x_polynomial {
            return false;
        }
        let degree = self.degree() as usize;
        let order = 1usize << degree - 1;
        let self_nonzero = self.clone();
        for q in prime_factors(order) {
            if x_polynomial.pow_mod(order / q, &self_nonzero) == *Self::one() {
                return false;
            }
        }
        true
    }

    /// Computes irreducible factors of square free polynomial using Berlekamp's algorithm.
    ///
    /// # Note
    /// If you know that your polynomial is square free (meaning all its irreducible factors are different), you can use [Self::square_free_irreducible_factors] instead.
    ///
    /// You can check whether a polynomial is square free using [Self::is_square_free].
    ///
    /// # Returns
    /// The set of irreducible factors and their respective power
    pub fn irreducible_factors(&self) -> HashMap<Self, usize> {
        // Thanks https://github.com/uranix/factormod/
        let mut polynomial_stack = Vec::from([(self.to_owned(), 1usize)]);
        let mut ret = HashMap::new();
        while let Some((polynomial, multiplier)) = polynomial_stack.pop() {
            let d = polynomial.double_factor();
            if d.is_one() {
                let factors = polynomial.square_free_irreducible_factors().unwrap();
                factors.iter().for_each(|factor| {
                    Self::incr_hashmap_count(&mut ret, factor, multiplier);
                });
                continue;
            }
            let poly_degree = polynomial.get().degree() as usize;
            if polynomial == d {
                // f' = 0, f only has even powers
                let mut g = BinaryPolynomial::zero();
                for i in (0..=poly_degree).step_by(2) {
                    if polynomial.get().polynomial.bit(i as u64) {
                        g.flip_bit(i >> 1);
                    }
                }
                let non_zero_g = NonZeroBinaryPolynomial::new(g).unwrap();
                polynomial_stack.push((non_zero_g, multiplier * 2));
                continue;
            }
            let add = NonZeroBinaryPolynomial::new(self.get().div_mod(&d).0)
                .unwrap()
                .square_free_irreducible_factors()
                .unwrap();
            add.iter().for_each(|p| {
                Self::incr_hashmap_count(&mut ret, p, multiplier);
            });
            polynomial_stack.push((d, multiplier));
        }
        ret
    }

    fn incr_hashmap_count(hashmap: &mut HashMap<Self, usize>, key: &Self, incr: usize) {
        if hashmap.contains_key(key) {
            hashmap.get_mut(key).unwrap().add_assign(incr);
        } else {
            hashmap.insert(key.to_owned(), incr);
        }
    }

    /// Computes irreducible factors of square free polynomial.
    ///
    /// A polynomial is square free if all its irreducible factors are differents.
    ///
    /// Compared to [Self::irreducible_factors], this function is faster but needs square free polynomial.
    ///
    /// You can check whether a polynomial is square free using [Self::is_square_free].
    ///
    /// # Returns
    /// The set of irreducible factors, or an error if the polynomial is not square free
    pub fn square_free_irreducible_factors(&self) -> Result<HashSet<Self>, BinaryPolynomialError> { // TODO better algorithm?
        if !self.is_square_free() {
            return Err(BinaryPolynomialError::NotSquareFreePolynomialError);
        }
        if self.is_one() {
            return Ok(HashSet::from_iter(vec![self.to_owned()]));
        }
        let mut b = BerkelampMatrix::from_poly(self);
        b.add_unit_matrix();
        let piv = b.rref();
        let hv = b.rref_nullspace(&piv);

        let mut hq = VecDeque::from(hv);
        assert!(!hq.is_empty());
        assert!(hq.front().unwrap().is_one());
        hq.pop_front();

        let mut in_set = HashSet::new();
        let mut out_set = HashSet::new();
        in_set.insert(self.to_owned());

        while let Some(polynomial) = hq.pop_front() {
            let h0 = NonZeroBinaryPolynomial::new(polynomial).unwrap();
            let mut h1 = h0.clone();
            h1.get_mut().flip_bit(0);
            out_set.clear();
            for p in &in_set {
                let d0 = p.gcd(&h0);
                let d1 = p.gcd(&h1);
                if !d0.is_one() {
                    out_set.insert(d0);
                }
                if !d1.is_one() {
                    out_set.insert(d1);
                }
            }
            std::mem::swap(&mut in_set, &mut out_set);
        }

        Ok(in_set)
    }

    fn double_factor(&self) -> Self {
        let derivative = self.get().derivative();
        let nonzero_derivative = match NonZeroBinaryPolynomial::new(derivative) {
            None => {
                return self.clone();
            }
            Some(d) => {d}
        };
        self.gcd(&nonzero_derivative)
    }

    /// Check if the polynomial is square free
    ///
    /// A polynomial is square free if all its irreducible factors are different
    ///
    /// # Returns
    /// `true` if the polynomial is square free, `false` otherwise
    pub fn is_square_free(&self) -> bool {
        self.double_factor().is_one()
    }
}

impl Deref for NonZeroBinaryPolynomial {
    type Target = BinaryPolynomial;

    fn deref(&self) -> &Self::Target {
        self.get()
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

/// Performs the multiplication of two `NonZeroBinaryPolynomial`
impl Mul for NonZeroBinaryPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.0 * rhs.0).unwrap()
    }
}

impl One for NonZeroBinaryPolynomial {
    fn one() -> Self {
        Self::new(BinaryPolynomial::one()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use crate::error::BinaryPolynomialError;
    use crate::{BinaryPolynomial, NonZeroBinaryPolynomial};
    use num_traits::{One, Zero};
    use std::collections::HashSet;
    use num_bigint::BigUint;

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

    #[test]
    fn test_square_free_irreducible_factors() {
        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2 + 1").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors();
        assert!(factorization.is_err());
        assert_eq!(factorization.unwrap_err(), BinaryPolynomialError::NotSquareFreePolynomialError); // (x + 1)^2

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x + 1").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors().unwrap();
        assert_eq!(factorization.len(), 1);
        assert_eq!(factorization.iter().map(|f| f.to_string()).collect::<Vec<String>>(), vec!["x + 1"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors().unwrap();
        assert_eq!(factorization.len(), 1);
        assert_eq!(factorization.iter().map(|f| f.to_string()).collect::<Vec<String>>(), vec!["x"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("1").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors().unwrap();
        assert_eq!(factorization.len(), 1);
        assert_eq!(factorization.iter().map(|f| f.to_string()).collect::<Vec<String>>(), vec!["1"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2 + x").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors().unwrap();
        assert_eq!(factorization.len(), 2);
        assert_eq!(factorization.iter().map(|p| p.to_string()).collect::<HashSet<String>>(), HashSet::from_iter(vec!["x + 1".to_string(), "x".to_string()]));

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^4 + x").unwrap()).unwrap();
        let factorization = p.square_free_irreducible_factors().unwrap();
        assert_eq!(factorization.len(), 3);
        assert_eq!(factorization.iter().map(|p| p.to_string()).collect::<HashSet<String>>(), HashSet::from_iter(vec!["x^2 + x + 1".to_string(), "x + 1".to_string(), "x".to_string()]));
    }

    #[test]
    fn test_irreducible_factors() {
        let p = NonZeroBinaryPolynomial::one();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<Vec<String>>(), vec!["1: 1"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<Vec<String>>(), vec!["x: 1"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x + 1").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<Vec<String>>(), vec!["x + 1: 1"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<Vec<String>>(), vec!["x: 2"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^3").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 1);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<Vec<String>>(), vec!["x: 3"]);

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2 + x").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 2);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<HashSet<String>>(), HashSet::from_iter(vec!["x: 1".to_string(), "x + 1: 1".to_string()]));

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^4 + x").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 3);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<HashSet<String>>(), HashSet::from_iter(vec!["x: 1".to_string(), "x + 1: 1".to_string(), "x^2 + x + 1: 1".to_string()]));

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^36 + x^34 + x^33 + x^4 + x^2 + x").unwrap()).unwrap();
        let factors = p.irreducible_factors();
        assert_eq!(factors.len(), 3);
        assert_eq!(factors.iter().map(|(p, n)| [p.to_string(), n.to_string()].join(": ")).collect::<HashSet<String>>(), HashSet::from_iter(vec!["x: 1".to_string(), "x + 1: 32".to_string(), "x^3 + x + 1: 1".to_string()]));
    }

    #[test]
    fn test_is_square_free() {
        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2 + x").unwrap()).unwrap();
        assert!(p.is_square_free());

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x + 1").unwrap()).unwrap();
        assert!(p.is_square_free());

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^4 + x").unwrap()).unwrap();
        assert!(p.is_square_free());

        let p = NonZeroBinaryPolynomial::new(BinaryPolynomial::try_from("x^2 + 1").unwrap()).unwrap();
        assert!(!p.is_square_free());
    }

    #[test]
    fn test_is_irreducible() {
        assert!(!NonZeroBinaryPolynomial::one().is_irreducible());
        assert!(!NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b11011usize))).unwrap().is_irreducible());
        assert!(NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b11111usize))).unwrap().is_irreducible());
        assert!(NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b10011usize))).unwrap().is_irreducible());
    }

    #[test]
    fn test_is_primitive() {
        assert!(!NonZeroBinaryPolynomial::one().is_primitive());
        assert!(!NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b11011usize))).unwrap().is_primitive());
        assert!(NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b11111usize))).unwrap().is_primitive());
        assert!(NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b10011usize))).unwrap().is_primitive());
    }
}
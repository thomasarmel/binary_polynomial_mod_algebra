//! # Binary univariate Polynomials
//! This crate provides a simple implementation of binary univariate polynomials and operations on them.

#![forbid(unsafe_code, unused_must_use)]
#![deny(
    missing_docs,
    unreachable_pub,
    unused_import_braces,
    unused_extern_crates,
    unused_qualifications
)]
#[doc = include_str!("../README.md")]

/// Errors in polynomial operations
pub mod error;
mod utils;
pub mod non_zero_polynomial;
mod berkelamp_matrix;

use crate::error::BinaryPolynomialError;
pub use crate::non_zero_polynomial::NonZeroBinaryPolynomial;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Pow, Zero};
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Rem};

/// Represents a binary univariate polynomial
#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd, Hash)]
pub struct BinaryPolynomial {
    polynomial: BigUint,
}

impl BinaryPolynomial {

    /// Returns the degree of the polynomial
    ///
    /// Degree of -1 is returned for the zero polynomial
    pub fn degree(&self) -> isize {
        self.polynomial.bits() as isize - 1
    }

    /// Performs the division of self-polynomial by the divisor polynomial
    ///
    /// # Arguments
    /// * `divisor` - The non-zero divisor polynomial
    ///
    /// # Returns
    /// A tuple of (quotient, remainder)
    pub fn div_mod(&self, divisor: &NonZeroBinaryPolynomial) -> (Self, Self) {
        let divisor = divisor.get();
        let mut q = BigUint::zero();
        let mut a = self.polynomial.clone();
        let bl = divisor.polynomial.bits() as i64;
        loop {
            let shift = (a.bits() as i64) - bl;
            if shift < 0 {
                return (Self { polynomial: q}, Self {polynomial: a});
            }
            q.set_bit(shift as u64, true);
            a ^= divisor.polynomial.clone() << shift;
        }
    }

    /// Performs the multiplication of self-polynomial by the other polynomial, modulo the modulus polynomial
    ///
    /// # Arguments
    /// * `other` - The other polynomial to multiply by
    /// * `modulus` - The non-zero modulus polynomial
    ///
    /// # Returns
    /// The result of the multiplication modulo the modulus polynomial, or an error if the degree of the multiplier or the other polynomial is greater or equal to the modulus polynomial
    pub fn mul_mod(&self, other: &Self, modulus: &NonZeroBinaryPolynomial) -> Result<Self, BinaryPolynomialError> {
        let modulus = modulus.get();
        let modulus_degree = modulus.degree();

        if self.degree() >= modulus_degree || other.degree() >= modulus_degree {
            return Err(BinaryPolynomialError::MultiplierDegreeGreaterOrEqualToModulusError);
        }
        let mut result = BigUint::zero();
        let mut a = self.polynomial.clone();
        let mut b = other.polynomial.clone();
        while !a.is_zero() && !b.is_zero() {
            if a.is_odd() {
                result ^= &b;
            }
            a >>= 1;
            b <<= 1;
            if (b.clone() >> modulus_degree).is_odd() {
                b ^= modulus.polynomial.clone();
            }
        }

        Ok(Self { polynomial: result })
    }

    /// Performs the exponentiation of the polynomial to the power of the exponent, modulo the modulus polynomial
    ///
    /// # Arguments
    /// * `exp` - The exponent to raise the polynomial to
    /// * `modulus` - The non-zero modulus polynomial
    ///
    /// # Returns
    /// The result of the exponentiation modulo the modulus polynomial
    pub fn pow_mod(&self, exp: usize, modulus: &NonZeroBinaryPolynomial) -> Self {
        let mut factor = self.clone() % modulus.clone();

        if self.is_zero() || self.is_one() || exp == 1 {
            return factor;
        }

        let mut result = Self::one();
        let mut exp = exp;
        while exp > 0 {
            if exp.is_odd() {
                result = result.mul_mod(&factor, modulus).unwrap();
            }
            factor = factor.mul_mod(&factor, modulus).unwrap();
            exp >>= 1;
        }
        result
    }

    /// Checks if the polynomial is congruent to the other polynomial modulo the modulus polynomial
    ///
    /// # Arguments
    /// * `other` - The other polynomial to check congruence with
    /// * `modulus` - The non-zero modulus polynomial
    ///
    /// # Returns
    /// `true` if the polynomials are congruent, `false` otherwise
    pub fn congruent_mod(&self, other: &Self, modulus: &NonZeroBinaryPolynomial) -> bool {
        (self.clone() + other.clone()) % modulus.clone() == Self::zero()
    }

    /// Returns the derivative of the polynomial
    pub fn derivative(&self) -> Self {
        let mut derivative = Self::zero();
        for i in (1..self.polynomial.bits()).step_by(2) {
            derivative.polynomial.set_bit(i-1, self.polynomial.bit(i));
        }
        derivative
    }

    fn flip_bit(&mut self, bit_pos: usize) {
        let bit_pos = bit_pos as u64;
        self.polynomial.set_bit(bit_pos, self.polynomial.bit(bit_pos) ^ true);
    }
}

/// Display implementation for `BinaryPolynomial`
///
/// The polynomial is displayed in the form of a sum of monomials, like: `x^3 + x^2 + 1`, or `0` for the zero polynomial
impl Display for BinaryPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            write!(f, "0")?;
            return Ok(());
        }
        let mut already_displayed = false;
        for i in (0..self.polynomial.bits()).rev() {
            let c = self.polynomial.bit(i);
            if c {
                if already_displayed {
                    write!(f, " + ")?;
                }
                already_displayed = true;
                if i == 0 {
                    write!(f, "1")?;
                } else if i == 1 {
                    write!(f, "x")?;
                }
                else {
                    write!(f, "x^{}", i)?;
                }
            }
        }
        Ok(())
    }
}

/// Conversion from BigUint to BinaryPolynomial
///
/// The polynomial is generated from the bits of the BigUint, with the least significant bit being lowest degree monoial
///
/// # Example
/// ```rust
/// use num_bigint::BigUint;
/// use binary_polynomial_mod_algebra::BinaryPolynomial;
/// let polynomial = BinaryPolynomial::from(BigUint::from(13u32)); // 0b1101
/// assert_eq!(polynomial.to_string(), "x^3 + x^2 + 1");
/// ```
impl From<BigUint> for BinaryPolynomial {
    fn from(polynomial: BigUint) -> Self {
        Self { polynomial }
    }
}

/// Conversion from Vec<bool> to BinaryPolynomial
///
/// The polynomial is generated from the bits of the Vec<bool>, with the first element being the highest degree monomial
///
/// # Example
/// ```rust
/// use binary_polynomial_mod_algebra::BinaryPolynomial;
/// let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
/// assert_eq!(polynomial.to_string(), "x^3 + x^2 + 1");
/// ```
impl From<Vec<bool>> for BinaryPolynomial {
    fn from(polynomial: Vec<bool>) -> Self {
        let mut result_polynomial = BigUint::zero();
        for (i, bit) in polynomial.iter().enumerate() {
            let bit_position = (polynomial.len() - i - 1) as u64;
            if *bit {
                result_polynomial.set_bit(bit_position, true);
            }
        }
        Self { polynomial: result_polynomial }
    }
}

impl TryFrom<&str> for BinaryPolynomial {
    type Error = BinaryPolynomialError;

    /// Parse polynomial from string
    ///
    /// # Arguments
    /// `value` - The string representation of the binary polynomial, in the form `x^4 + x^3 + x + 1`
    ///
    /// # Returns
    /// A `BinaryPolynomial`, or an error if the parsing failed
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let polynomial_filtered = value.chars().filter(|c| !c.is_whitespace()).collect::<String>();
        if polynomial_filtered.is_empty() {
            return Ok(Self::zero());
        }
        let mut polynomial = Self::zero();
        for monomial_string in polynomial_filtered.split('+') {
            if monomial_string == "0" {
                continue;
            }
            if monomial_string == "1" {
                polynomial.flip_bit(0);
                continue;
            }
            if monomial_string == "x" {
                polynomial.flip_bit(1);
                continue;
            }
            let monomial_splitted = monomial_string.split('^').collect::<Vec<&str>>();
            if monomial_splitted.len() != 2 || monomial_splitted[0] != "x" {
                return Err(BinaryPolynomialError::ParsingPolynomialError)
            }
            let monomial_pos = monomial_splitted[1].parse::<usize>().map_err(|_| BinaryPolynomialError::ParsingPolynomialError)?;
            polynomial.flip_bit(monomial_pos);
        }
        Ok(polynomial)
    }
}

/// Conversion from BinaryPolynomial to BigUint
///
/// The BigUint is generated from the bits of the BinaryPolynomial, with the least significant bit being lowest degree monomial
impl Into<BigUint> for BinaryPolynomial {
    fn into(self) -> BigUint {
        self.polynomial
    }
}

/// Conversion from BinaryPolynomial to Vec<bool>
///
/// The Vec<bool> is generated from the bits of the BinaryPolynomial, with the first element being the highest degree monomial
impl Into<Vec<bool>> for BinaryPolynomial {
    fn into(self) -> Vec<bool> {
        let mut result = Vec::new();
        for i in (0..self.polynomial.bits()).rev() {
            result.push(self.polynomial.bit(i));
        }
        result
    }
}

/// Performs the addition of two BinaryPolynomials
///
/// This is equivalent to the XOR operation on the bits of the polynomials
impl Add for BinaryPolynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            polynomial: &self.polynomial ^ &rhs.polynomial,
        }
    }
}

impl Zero for BinaryPolynomial {
    fn zero() -> Self {
        Self {
            polynomial: BigUint::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.polynomial.is_zero()
    }
}

/// Performs the multiplication of two `BinaryPolynomial`
impl Mul<Self> for BinaryPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut a = self.polynomial.clone();
        let mut b = rhs.polynomial.clone();
        let mut result = BigUint::zero();
        while !a.is_zero() && !b.is_zero() {
            if a.is_odd() {
                result ^= &b;
            }
            a >>= 1;
            b <<= 1;
        }
        Self {
            polynomial: result,
        }
    }
}

impl One for BinaryPolynomial {
    fn one() -> Self {
        Self {
            polynomial: BigUint::one(),
        }
    }
}

/// Performs the remainder operation of a binary polynomial by a non-zero binary polynomial
impl Rem<NonZeroBinaryPolynomial> for BinaryPolynomial {
    type Output = Self;

    fn rem(self, rhs: NonZeroBinaryPolynomial) -> Self::Output {
        let rhs = rhs.get();
        let bl = rhs.polynomial.bits() as i64;
        let mut a = self.polynomial.clone();
        loop {
            let shift = (a.bits() as i64) - bl;
            if shift < 0 {
                return Self { polynomial: a };
            }
            a ^= rhs.polynomial.clone() << shift;
        }
    }
}

/// Performs the division of self-polynomial by the divisor polynomial
impl Div<NonZeroBinaryPolynomial> for BinaryPolynomial {
    type Output = Self;

    fn div(self, rhs: NonZeroBinaryPolynomial) -> Self::Output {
        self.div_mod(&rhs).0
    }
}

/// Performs the exponentiation of a binary polynomial to an usize exponent
impl Pow<usize> for BinaryPolynomial {
    type Output = Self;

    fn pow(self, exp: usize) -> Self::Output {
        if self.is_zero() || self.is_one() || exp == 1 {
            return self;
        }
        if exp == 0 {
            return Self::one();
        }
        let mut factor = self.clone();
        let mut result = Self::one();
        let mut exp = exp;
        while exp > 0 {
            if exp.is_odd() {
                result = result * factor.clone();
            }
            factor = factor.clone() * factor;
            exp >>= 1;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display() {
        let zero_polynomial = BinaryPolynomial::zero();
        assert_eq!(zero_polynomial.to_string(), "0");

        let one_polynomial = BinaryPolynomial::one();
        assert_eq!(one_polynomial.to_string(), "1");

        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        assert_eq!(polynomial.to_string(), "x^2 + 1");

        let polynomial = BinaryPolynomial::from(vec![true, false, true, false]);
        assert_eq!(polynomial.to_string(), "x^3 + x");
    }

    #[test]
    fn test_is_zero() {
        let zero_polynomial = BinaryPolynomial::zero();
        assert!(zero_polynomial.is_zero());

        let polynomial = BinaryPolynomial::from(vec![false, false, false, false]);
        assert!(polynomial.is_zero());

        let one_polynomial = BinaryPolynomial::one();
        assert!(!one_polynomial.is_zero());

        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        assert!(!polynomial.is_zero());
    }

    #[test]
    fn test_from_biguint() {
        let polynomial = BinaryPolynomial::from(BigUint::from(0b1010u32));
        assert_eq!(polynomial.to_string(), "x^3 + x");
        let biguint: BigUint = polynomial.into();
        assert_eq!(biguint, BigUint::from(0b1010u32));

        let polynomial = BinaryPolynomial::from(BigUint::from(0b1101u32));
        assert_eq!(polynomial.to_string(), "x^3 + x^2 + 1");
        let biguint: BigUint = polynomial.into();
        assert_eq!(biguint, BigUint::from(0b1101u32));
    }

    #[test]
    fn test_from_vec_bool() {
        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        assert_eq!(polynomial.to_string(), "x^2 + 1");
        let vec_bool: Vec<bool> = polynomial.into();
        assert_eq!(vec_bool, vec![true, false, true]);

        let polynomial = BinaryPolynomial::from(vec![true, false, true, false]);
        assert_eq!(polynomial.to_string(), "x^3 + x");
        let vec_bool: Vec<bool> = polynomial.into();
        assert_eq!(vec_bool, vec![true, false, true, false]);
    }

    #[test]
    fn test_from_str() {
        let polynomial = BinaryPolynomial::try_from("").unwrap();
        assert!(polynomial.is_zero());
        let polynomial = BinaryPolynomial::try_from("  ").unwrap();
        assert!(polynomial.is_zero());
        let polynomial = BinaryPolynomial::try_from(" 0").unwrap();
        assert!(polynomial.is_zero());
        let polynomial = BinaryPolynomial::try_from("1 ").unwrap();
        assert!(polynomial.is_one());
        let polynomial = BinaryPolynomial::try_from("x + x^31 + 1").unwrap();
        assert_eq!(polynomial.to_string(), "x^31 + x + 1");
        let polynomial = BinaryPolynomial::try_from("x^31 + x^ 32 + x + x^31 + 1 + 1").unwrap();
        assert_eq!(polynomial.to_string(), "x^32 + x");
        let polynomial = BinaryPolynomial::try_from("x^3 + a^2 + 1");
        assert!(polynomial.is_err());
        assert_eq!(polynomial.unwrap_err(), BinaryPolynomialError::ParsingPolynomialError);
        let polynomial = BinaryPolynomial::try_from("x^3 + x^a + 1");
        assert!(polynomial.is_err());
        assert_eq!(polynomial.unwrap_err(), BinaryPolynomialError::ParsingPolynomialError);
    }

    #[test]
    fn test_add() {
        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true, true]);
        let result = polynomial + polynomial2;
        assert_eq!(result.to_string(), "x^3 + x^2 + x");
        let vec_bool: Vec<bool> = result.into();
        assert_eq!(vec_bool, vec![true, true, true, false]);
    }

    #[test]
    fn test_mul() {
        let zero_polynomial = BinaryPolynomial::zero();
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true, true]);

        let result = zero_polynomial * polynomial.clone();
        assert!(result.is_zero());

        let result = one_polynomial * polynomial.clone();
        assert_eq!(result, polynomial);

        let result = polynomial.clone() * polynomial2.clone();
        assert_eq!(result.to_string(), "x^5 + x^2 + x + 1");
    }

    #[test]
    fn test_modulo() {
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let modulo = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, true])).unwrap();

        let result = one_polynomial.clone() % modulo.clone();
        assert_eq!(result, one_polynomial);

        let result = polynomial % modulo;
        assert_eq!(result.to_string(), "x + 1");
    }

    #[test]
    fn test_div_mod() {
        let zero_polynomial = BinaryPolynomial::zero();
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);

        let result = polynomial.div_mod(&NonZeroBinaryPolynomial::new(one_polynomial).unwrap());
        assert_eq!(result, (polynomial.clone(), zero_polynomial));

        let result = polynomial.div_mod(&NonZeroBinaryPolynomial::new(polynomial2).unwrap());
        assert_eq!(result.0.to_string(), "x^2 + x");
        assert_eq!(result.1.to_string(), "x + 1");
    }

    #[test]
    fn test_div() {
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);

        let result = polynomial.clone() / NonZeroBinaryPolynomial::new(one_polynomial).unwrap();
        assert_eq!(result, polynomial.clone());

        let result = polynomial / NonZeroBinaryPolynomial::new(polynomial2).unwrap();
        assert_eq!(result.to_string(), "x^2 + x");
    }

    #[test]
    fn test_degree() {
        let zero_polynomial = BinaryPolynomial::zero();
        assert_eq!(zero_polynomial.degree(), -1);

        let one_polynomial = BinaryPolynomial::one();
        assert_eq!(one_polynomial.degree(), 0);

        let x_polynomial = BinaryPolynomial::from(vec![true, false]);
        assert_eq!(x_polynomial.degree(), 1);

        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        assert_eq!(polynomial.degree(), 4);
    }

    #[test]
    fn test_mul_mod() {
        let zero_polynomial = BinaryPolynomial::zero();
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let polynomial3 = BinaryPolynomial::from(vec![true, true, false, true]);
        let modulo = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, false, true])).unwrap();

        let result = zero_polynomial.mul_mod(&polynomial2, &modulo).unwrap();
        assert_eq!(result.to_string(), "0");

        let result = one_polynomial.mul_mod(&polynomial2, &modulo).unwrap();
        assert_eq!(result, polynomial2);

        let result = polynomial.mul_mod(&polynomial2, &modulo);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::MultiplierDegreeGreaterOrEqualToModulusError);

        let result = polynomial2.mul_mod(&polynomial, &modulo);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::MultiplierDegreeGreaterOrEqualToModulusError);

        let result = polynomial2.mul_mod(&polynomial3, &modulo).unwrap();
        assert_eq!(result.to_string(), "x^3 + x");

        let result = polynomial3.mul_mod(&polynomial2, &modulo).unwrap();
        assert_eq!(result.to_string(), "x^3 + x");

        let result = polynomial3.mul_mod(&polynomial3, &modulo).unwrap();
        assert_eq!(result.to_string(), "x^2");
    }

    #[test]
    fn test_pow() {
        let zero_polynomial = BinaryPolynomial::zero();
        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);

        let result = zero_polynomial.pow(3);
        assert!(result.is_zero());

        let result = polynomial.clone().pow(0);
        assert!(result.is_one());

        let result = polynomial.clone().pow(1);
        assert_eq!(result, polynomial);

        let result = polynomial.clone().pow(2);
        assert_eq!(result.to_string(), "x^6 + x^4 + 1");

        let result = polynomial.clone().pow(3);
        assert_eq!(result.to_string(), "x^9 + x^8 + x^7 + x^4 + x^3 + x^2 + 1");
    }

    #[test]
    fn test_pow_mod() {
        let zero_polynomial = BinaryPolynomial::zero();
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
        let modulo = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, false, true])).unwrap();

        let result = zero_polynomial.pow_mod(3, &modulo);
        assert!(result.is_zero());

        let result = one_polynomial.pow_mod(3, &modulo);
        assert_eq!(result, one_polynomial);

        let result = polynomial.pow_mod(0, &modulo);
        assert_eq!(result, one_polynomial);

        let result = polynomial.pow_mod(1, &modulo);
        assert_eq!(result, polynomial);

        let result = polynomial.pow_mod(2, &modulo);
        assert_eq!(result.to_string(), "x^2");

        let result = polynomial.pow_mod(3, &modulo);
        assert_eq!(result.to_string(), "x^2 + x + 1");
    }

    #[test]
    fn test_congruent_mod() {
        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let modulo = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(vec![true, false, false, false, true])).unwrap();

        let result = polynomial.congruent_mod(&polynomial2, &modulo);
        assert!(!result);

        let result = polynomial.congruent_mod(&polynomial, &modulo);
        assert!(result);
    }

    #[test]
    fn test_derivative() {
        assert_eq!(BinaryPolynomial::zero().derivative(), BinaryPolynomial::zero());
        assert_eq!(BinaryPolynomial::one().derivative(), BinaryPolynomial::zero());
        assert_eq!(BinaryPolynomial::from(BigUint::from(0b10usize)).derivative(), BinaryPolynomial::one());
        assert_eq!(BinaryPolynomial::from(BigUint::from(0b100usize)).derivative(), BinaryPolynomial::zero());
        assert_eq!(BinaryPolynomial::from(BigUint::from(0b110usize)).derivative(), BinaryPolynomial::one());
    }
}

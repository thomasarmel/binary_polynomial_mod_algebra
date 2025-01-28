pub mod error;

use std::fmt::Display;
use std::ops::{Add, Mul, Rem};
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, Pow, Zero};
use crate::error::BinaryPolynomialError;

#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd)]
pub struct BinaryPolynomial {
    polynomial: BigUint,
}

impl BinaryPolynomial {

    pub fn degree(&self) -> usize {
        if self.polynomial.is_zero() {
            return 0;
        }
        self.polynomial.bits() as usize - 1
    }

    /// Returns (quotient, remainder)
    pub fn div_mod(&self, divisor: &Self) -> Result<(Self, Self), BinaryPolynomialError> {
        if divisor.is_zero() {
            return Err(BinaryPolynomialError::DivideByZeroError);
        }
        let mut q = BigUint::zero();
        let mut a = self.polynomial.clone();
        let bl = divisor.polynomial.bits() as i64;
        loop {
            let shift = (a.bits() as i64) - bl;
            if shift < 0 {
                return Ok(
                    (Self { polynomial: q}, Self {polynomial: a})
                );
            }
            q ^= BigUint::one() << shift;
            a ^= divisor.polynomial.clone() << shift;
        }
    }

    pub fn mul_mod(&self, other: &Self, modulus: &Self) -> Result<Self, BinaryPolynomialError> {
        if modulus.is_zero() {
            return Err(BinaryPolynomialError::NullModulusError);
        }
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

    pub fn gcd(&self, other: &Self) -> Result<Self, BinaryPolynomialError> {
        if self.is_zero() || other.is_zero() {
            return Err(BinaryPolynomialError::NullPolynomialCommonDivisorError);
        }

        let mut a = self.clone();
        let mut b = other.clone();

        while !b.polynomial.is_zero() {
            (a, b) = (b.clone(), a.clone() % b.clone());
        }
        Ok(a)
    }

    /// Returns (d, x, y) where d is the Greatest Common Divisor of polynomials a and b, x, y are polynomials that satisfy: p_mul(a,x) ^ p_mul(b,y) = d
    pub fn egcd(&self, other: &Self) -> Result<(Self, Self, Self), BinaryPolynomialError> {
        if self.is_zero() || other.is_zero() {
            return Err(BinaryPolynomialError::NullPolynomialCommonDivisorError);
        }

        let mut a = (self.clone(), BinaryPolynomial::one(), BinaryPolynomial::zero());
        let mut b = (other.clone(), BinaryPolynomial::zero(), BinaryPolynomial::one());

        loop {
            let (q, r) = a.0.div_mod(&b.0)?;
            let q: BinaryPolynomial = q.into();
            let r:BinaryPolynomial = r.into();
            if r.is_zero() {
                return Ok(b);
            }
            (a, b) = (b.clone(), (r, BinaryPolynomial::from(a.1.polynomial ^ (q.clone() * b.1.clone()).polynomial), BinaryPolynomial::from(a.2.polynomial ^ (q.clone() * b.2.clone()).polynomial)));
        }
    }

    pub fn inv_mod(&self, modulus: &Self) -> Result<Self, BinaryPolynomialError> {
        if modulus.is_zero() {
            return Err(BinaryPolynomialError::NullModulusError);
        }
        if self.is_zero() {
            return Err(BinaryPolynomialError::NonInvertiblePolynomialError);
        }
        let (d, x, _) = Self::egcd(&self, modulus)?;
        if !d.is_one() {
            return Err(BinaryPolynomialError::NonInvertiblePolynomialError);
        }
        Ok(x)
    }

    pub fn pow_mod(&self, exp: usize, modulus: &Self) -> Result<Self, BinaryPolynomialError> {
        if modulus.is_zero() {
            return Err(BinaryPolynomialError::NullModulusError);
        }
        if self.is_zero() || self.is_one() || exp == 1 {
            return Ok(self.clone() % modulus.clone());
        }

        let mut factor = self.clone() % modulus.clone();
        let mut result = Self::one();
        let mut exp = exp;
        while exp > 0 {
            if exp.is_odd() {
                result = result.mul_mod(&factor, modulus)?;
            }
            factor = factor.mul_mod(&factor, modulus)?;
            exp >>= 1;
        }
        Ok(result)
    }

    pub fn congruent_mod(&self, other: &Self, modulus: &Self) -> bool {
        if modulus.is_zero() {
            panic!("Modulus cannot be zero");
        }
        (self.clone() + other.clone()) % modulus.clone() == Self::zero()
    }

    pub fn coprime(&self, other: &Self) -> bool {
        self.gcd(other).unwrap().is_one()
    }
}

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

impl From<BigUint> for BinaryPolynomial {
    fn from(polynomial: BigUint) -> Self {
        Self { polynomial }
    }
}

/// Biggest factor is first
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

impl Into<BigUint> for BinaryPolynomial {
    fn into(self) -> BigUint {
        self.polynomial
    }
}

/// Biggest factor is first
impl Into<Vec<bool>> for BinaryPolynomial {
    fn into(self) -> Vec<bool> {
        let mut result = Vec::new();
        for i in (0..self.polynomial.bits()).rev() {
            result.push(self.polynomial.bit(i));
        }
        result
    }
}

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

impl Rem for BinaryPolynomial {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            panic!("Modulus cannot be zero");
        }
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
            factor = factor.clone() * factor.clone();
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
        let modulo = BinaryPolynomial::from(vec![true, false, true]);

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

        let result = polynomial.div_mod(&zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::DivideByZeroError);

        let result = polynomial.div_mod(&one_polynomial).unwrap();
        assert_eq!(result, (polynomial.clone(), zero_polynomial));

        let result = polynomial.div_mod(&polynomial2).unwrap();
        assert_eq!(result.0.to_string(), "x^2 + x");
        assert_eq!(result.1.to_string(), "x + 1");
    }

    #[test]
    fn test_degree() {
        let zero_polynomial = BinaryPolynomial::zero();
        assert_eq!(zero_polynomial.degree(), 0);

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
        let modulo = BinaryPolynomial::from(vec![true, false, false, false, true]);

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

        let result = polynomial2.mul_mod(&polynomial3, &zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullModulusError);

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
    fn test_gcd() {
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let zero_polynomial = BinaryPolynomial::zero();

        let result = polynomial.gcd(&polynomial2).unwrap();
        assert_eq!(result.to_string(), "x + 1");

        let result = polynomial2.gcd(&polynomial).unwrap();
        assert_eq!(result.to_string(), "x + 1");

        let result = polynomial.gcd(&polynomial).unwrap();
        assert_eq!(result, polynomial);

        let result = polynomial.gcd(&zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullPolynomialCommonDivisorError);

        let result = zero_polynomial.gcd(&polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullPolynomialCommonDivisorError);
    }

    #[test]
    fn test_egcd() {
        let polynomial = BinaryPolynomial::from(vec![true, true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let zero_polynomial = BinaryPolynomial::zero();

        let result = polynomial.egcd(&zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullPolynomialCommonDivisorError);

        let result = zero_polynomial.egcd(&polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullPolynomialCommonDivisorError);

        let result = polynomial.egcd(&polynomial2).unwrap();
        assert_eq!(result.0.to_string(), "x + 1");
        assert_eq!(result.1.to_string(), "1");
        assert_eq!(result.2.to_string(), "x^2 + x");

        let result = polynomial2.egcd(&polynomial).unwrap();
        assert_eq!(result.0.to_string(), "x + 1");
        assert_eq!(result.1.to_string(), "x^2 + x");
        assert_eq!(result.2.to_string(), "1");

        let result = polynomial.egcd(&polynomial).unwrap();
        assert_eq!(result.0.to_string(), "x^4 + x^3 + x^2 + 1");
        assert_eq!(result.1.to_string(), "0");
        assert_eq!(result.2.to_string(), "1");
    }

    #[test]
    fn test_inv_mod() {
        let zero_polynomial = BinaryPolynomial::zero();
        let modulus_polynomial = BinaryPolynomial::from(vec![true, false, false, false, true]);
        let polynomial = BinaryPolynomial::from(vec![true, false, true, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let polynomial3 = BinaryPolynomial::from(vec![true, false, false, false]);

        let result = polynomial.inv_mod(&zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullModulusError);

        let result = zero_polynomial.inv_mod(&modulus_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NonInvertiblePolynomialError);

        let result = polynomial.inv_mod(&modulus_polynomial).unwrap();
        assert_eq!(result.to_string(), "x^3 + x + 1");

        let result = polynomial2.inv_mod(&modulus_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NonInvertiblePolynomialError);

        let result = polynomial3.inv_mod(&modulus_polynomial).unwrap();
        assert_eq!(result.to_string(), "x");
    }

    #[test]
    fn test_pow_mod() {
        let zero_polynomial = BinaryPolynomial::zero();
        let one_polynomial = BinaryPolynomial::one();
        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
        let modulo = BinaryPolynomial::from(vec![true, false, false, false, true]);

        let result = zero_polynomial.pow_mod(3, &modulo).unwrap();
        assert!(result.is_zero());

        let result = one_polynomial.pow_mod(3, &modulo).unwrap();
        assert_eq!(result, one_polynomial);

        let result = polynomial.pow_mod(0, &modulo).unwrap();
        assert_eq!(result, one_polynomial);

        let result = polynomial.pow_mod(3, &zero_polynomial);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), BinaryPolynomialError::NullModulusError);

        let result = polynomial.pow_mod(1, &modulo).unwrap();
        assert_eq!(result, polynomial);

        let result = polynomial.pow_mod(2, &modulo).unwrap();
        assert_eq!(result.to_string(), "x^2");

        let result = polynomial.pow_mod(3, &modulo).unwrap();
        assert_eq!(result.to_string(), "x^2 + x + 1");
    }

    #[test]
    fn test_congruent_mod() {
        let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, true]);
        let modulo = BinaryPolynomial::from(vec![true, false, false, false, true]);

        let result = polynomial.congruent_mod(&polynomial2, &modulo);
        assert!(!result);

        let result = polynomial.congruent_mod(&polynomial, &modulo);
        assert!(result);
    }

    #[test]
    fn test_coprime() {
        let polynomial = BinaryPolynomial::from(vec![false, true, false, true]);
        let polynomial2 = BinaryPolynomial::from(vec![true, false, false, true]);
        let polynomial3 = BinaryPolynomial::from(vec![true, true, false, true]);

        let result = polynomial.coprime(&polynomial2);
        assert!(!result);

        let result = polynomial.coprime(&polynomial);
        assert!(!result);

        let result = polynomial2.coprime(&polynomial3);
        assert!(result);
    }
}

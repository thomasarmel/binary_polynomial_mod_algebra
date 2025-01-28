#[derive(thiserror::Error, Debug, PartialEq, Eq, Copy, Clone)]
pub enum BinaryPolynomialError {
    #[error("Polynomial is null")]
    NullPolynomialError,
    #[error("Division by zero")]
    DivideByZeroError,
    #[error("Modulus cannot be zero")]
    NullModulusError,
    #[error("Multiplier degree cannot be greater or equal to modulus")]
    MultiplierDegreeGreaterOrEqualToModulusError,
    #[error("Cannot find common divisor of null polynomials")]
    NullPolynomialCommonDivisorError,
    #[error("Non-invertible polynomial")]
    NonInvertiblePolynomialError,
}
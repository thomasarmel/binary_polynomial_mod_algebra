//! # Error module

/// Error type for binary polynomial operations
#[derive(thiserror::Error, Debug, PartialEq, Eq, Copy, Clone)]
pub enum BinaryPolynomialError {
    /// Error when the degree of the polynomial is greater than the modulus on polynomial modular multiplication
    #[error("Multiplier degree cannot be greater or equal to modulus")]
    MultiplierDegreeGreaterOrEqualToModulusError,
    /// Error when the given polynomial is not square free, ie it can be expressed as A\*A\*B, with A, B polynomials
    #[error("The polynomial is not square free")]
    NotSquareFreePolynomialError,
    /// Error when parsing polynomial from string
    #[error("The string cannot be parsed as a polynomial")]
    ParsingPolynomialError,
}
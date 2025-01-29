//! # Error module

/// Error type for binary polynomial operations
#[derive(thiserror::Error, Debug, PartialEq, Eq, Copy, Clone)]
pub enum BinaryPolynomialError {
    /// Error when the degree of the polynomial is greater than the modulus on polynomial modular multiplication
    #[error("Multiplier degree cannot be greater or equal to modulus")]
    MultiplierDegreeGreaterOrEqualToModulusError,
}
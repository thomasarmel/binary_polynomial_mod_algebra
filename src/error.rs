#[derive(thiserror::Error, Debug, PartialEq, Eq, Copy, Clone)]
pub enum BinaryPolynomialError {
    #[error("Multiplier degree cannot be greater or equal to modulus")]
    MultiplierDegreeGreaterOrEqualToModulusError,
    #[error("Non-invertible polynomial")]
    NonInvertiblePolynomialError,
}
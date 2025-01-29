# Binary univariate polynomial algebra

*This repo is highly inspired from https://gist.github.com/mildsunrise/e21ae2b1649532813f2594932f9e9371*

---

## Installation

On your project's directory, run the following command:

```bash
cargo add binary_polynomial_mod_algebra
```

## Usage

### Polynomial creation

You can create a `BinaryPolynomial` from `BigUint`, where the least significant bit is the lowest degree coefficient.

```rust
use num_bigint::BigUint;
use binary_polynomial_mod_algebra::BinaryPolynomial;

let polynomial = BinaryPolynomial::from(BigUint::from(13u32)); // or 0b1101u32
assert_eq!(polynomial.to_string(), "x^3 + x^2 + 1");
```

You can also create a `BinaryPolynomial` from a `Vec<bool>`, where the first element is the highest degree coefficient.

```rust
use binary_polynomial_mod_algebra::BinaryPolynomial;

let polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
assert_eq!(polynomial.to_string(), "x^3 + x^2 + 1");
```

Finally, you can create the null or unitary BinaryPolynomial:

```rust
use binary_polynomial_mod_algebra::BinaryPolynomial;
use num_traits::identities::One;
use num_traits::Zero;

let polynomial = BinaryPolynomial::zero();
assert_eq!(polynomial.to_string(), "0");

let polynomial = BinaryPolynomial::one();
assert_eq!(polynomial.to_string(), "1");
```

### Non-zero polynomial

Some functions require a `NonZeroBinaryPolynomial` argument, you can create it from a `BinaryPolynomial`:

```rust
use binary_polynomial_mod_algebra::BinaryPolynomial;
use binary_polynomial_mod_algebra::NonZeroBinaryPolynomial;
use num_traits::identities::Zero;

let zero_polynomial = BinaryPolynomial::zero();
let try_non_zero_polynomial = NonZeroBinaryPolynomial::new(zero_polynomial);
assert!(try_non_zero_polynomial.is_none());

let non_zero_polynomial = BinaryPolynomial::from(vec![true, true, false, true]);
let non_zero_polynomial = NonZeroBinaryPolynomial::new(non_zero_polynomial).unwrap();
assert!(non_zero_polynomial.get().to_string() == "x^3 + x^2 + 1");
```

#### Example

Polynomial inversion modulo a non-zero polynomial:

```rust
use binary_polynomial_mod_algebra::BinaryPolynomial;
use binary_polynomial_mod_algebra::NonZeroBinaryPolynomial;

let polynomial = NonZeroBinaryPolynomial::new(
        BinaryPolynomial::from(vec![true, false, true, true]) // x^3 + x + 1
    ).unwrap();
let modulo = NonZeroBinaryPolynomial::new(
        BinaryPolynomial::from(vec![true, false, false, false, true]) // x^4 + 1
    ).unwrap();
let inverse = polynomial.inv_mod(&modulo);
assert!(inverse.is_some()); // Inverse exists
assert_eq!(inverse.unwrap().to_string(), "x^3 + x + 1");
```
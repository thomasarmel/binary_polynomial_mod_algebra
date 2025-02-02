use std::collections::HashSet;
use num_integer::Integer;

pub(super) fn prime_factors(num: usize) -> HashSet<usize> {
    let mut factors = HashSet::new();
    if num == 0 {
        return factors;
    }
    let mut n = num;
    if n.is_even() {
        factors.insert(2);
        n >>= n.trailing_zeros() as usize;
    }

    let mut i = 3usize;
    while i*i <= n {
        while n % i == 0 {
            factors.insert(i);
            n /= i;
        }
        i += 2;
    }
    if n > 2 {
        factors.insert(n);
    }
    factors
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use crate::utils::prime_factors;

    #[test]
    fn test_prime_factors() {
        assert_eq!(prime_factors(0), HashSet::from([]));
        assert_eq!(prime_factors(1), HashSet::from([]));
        assert_eq!(prime_factors(7), HashSet::from([7]));
        assert_eq!(prime_factors(10), HashSet::from([5, 2]));
    }
}
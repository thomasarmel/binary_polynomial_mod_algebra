use std::fmt::{Display, Formatter};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use crate::{BinaryPolynomial, NonZeroBinaryPolynomial};

#[derive(Debug, Clone)]
pub(super) struct BerkelampMatrix {
    data: Vec<Vec<bool>>,
}

impl BerkelampMatrix {
    pub(super) fn from_poly(poly: &NonZeroBinaryPolynomial) -> Self {
        let poly_degree = poly.get().degree() as usize;
        let mut b = Self {
            data: vec![vec![false; poly_degree]; poly_degree],
        };
        let mut q = BinaryPolynomial::one();
        let x2 = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b100usize))).unwrap();
        for i in 0..poly_degree {
            for j in 0..poly_degree {
                b.data[j][i] ^= q.polynomial.bit(j as u64);
            }
            q = (q * x2.clone().get_owned()) % poly.clone();
        }
        b
    }

    pub(super) fn add_unit_matrix(&mut self) {
        for i in 0..self.data.len() {
            self.data[i][i] ^= true;
        }
    }

    pub(super) fn rref(&mut self) -> Vec<usize> {
        let matrix_len = self.data.len();
        let mut piv = Vec::<usize>::with_capacity(matrix_len);
        for j in 0..matrix_len {
            let mut pivrow = None;
            for i in piv.len()..matrix_len {
                if self.data[i][j] {
                    pivrow = Some(i);
                    break;
                }
            }
            if pivrow.is_none() {
                continue;
            }
            let mut pivrow = pivrow.unwrap();
            let i = piv.len();
            if i != pivrow {
                self.data.swap(i, pivrow);
            }
            piv.push(j);
            pivrow = i;
            for i in 0..matrix_len {
                if i != pivrow {
                    if self.data[i][j] {
                        for pos in 0..matrix_len {
                            self.data[i][pos] ^= self.data[pivrow][pos];
                        }
                    }
                }
            }
        }
        piv
    }

    pub(super) fn rref_nullspace(&self, piv: &[usize]) -> Vec<BinaryPolynomial> {
        let piv_len = piv.len();
        let matrix_size = self.data.len();
        let mut ret: Vec<BinaryPolynomial> = Vec::with_capacity(piv_len);
        let mut k = 0usize;
        for j in 0..matrix_size {
            if k < piv_len && j == piv[k] {
                k+=1;
                continue;
            }
            let mut p = BinaryPolynomial::zero();
            p.flip_bit(j);
            for i in 0..piv_len {
                if self.data[i][j] {
                    p.flip_bit(piv[i]);
                }
            }
            ret.push(p);
        }
        ret
    }
}

impl Display for BerkelampMatrix {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.data.len() {
            for j in 0..self.data[i].len() {
                write!(f, "{}", if self.data[i][j] {'1'}else{'0'})?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use crate::{BinaryPolynomial, NonZeroBinaryPolynomial};
    use crate::berkelamp_matrix::BerkelampMatrix;

    #[test]
    fn test_from_poly() {
        let poly = NonZeroBinaryPolynomial::new(BinaryPolynomial::from(BigUint::from(0b10110usize))).unwrap();
        let m = BerkelampMatrix::from_poly(&poly);
        assert_eq!(m.data, [[true, false, false, false], [false, false, true, true], [false, true, true, true], [false, false, false, true]]);
    }
}
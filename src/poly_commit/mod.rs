pub mod errors;
pub mod transcript;
pub mod gens;

use curve25519_dalek::ristretto::RistrettoPoint;
use crate::poly_commit::errors::ProofResult;
use curve25519_dalek::scalar::Scalar;
use crate::types::GroupElt;
use polynomials::Polynomial;

#[cfg(feature = "std")]
use num_bigint::*;

#[cfg(feature = "std")]
use rand::{CryptoRng, RngCore};

pub struct CommitmentKeyPC {
    pp: gens::ProofGens,
    H: super::types::GroupElt,
}

pub trait PolyCommitment {
    fn commit(&self, p: Polynomial<Scalar>, deg: u32, r: Scalar) -> ProofResult<GroupElt>;
    fn open<R: RngCore + CryptoRng>(&self, p_poly: Polynomial<Scalar>, commitment: GroupElt, deg: u32, z: Scalar, w: Scalar, rng: &mut R) -> ProofResult<()>;
}

pub struct PolyCommitmentScheme {
    ck: CommitmentKeyPC,
}

impl PolyCommitment for PolyCommitmentScheme {
    fn commit(&self, p: Polynomial<Scalar>, deg: u32, r: Scalar) -> ProofResult<GroupElt> {
        assert!(p.degree() == deg as usize);
        let commit = self.ck.pp.commit(p.into(), r);
        commit
    }

    #[cfg(feature = "std")]
    fn open<R: RngCore + CryptoRng>(&self, p: Polynomial<Scalar>, commitment: GroupElt, deg: u32, z: Scalar, w: Scalar, rng: &mut R) -> ProofResult<()> {
        let p_eval = p.eval(z);
        // sample random poly with p_bar(z) = 0
        let p_bar = {
            let mut partial_poly: Polynomial<Scalar> = Polynomial::new();
            for i in 0..((deg - 1) as usize) {
                partial_poly.push(Scalar::random(rng));
            }

            // z_root_polly = (x - z)
            let mut z_root_poly: Polynomial<Scalar> = Polynomial::new();
            z_root_poly.push(z);
            z_root_poly.push(Scalar::from(1u32));

            partial_poly * z_root_poly
        };

        let w_bar = Scalar::random(rng);
        // compute hiding commitment: \hat{C} = CM.Commit(ck, coeff_p_bar, w_bar)
        let hiding_commitment_C_bar = {
            let coeff_p_bar: Vec<Scalar> = p_bar.clone().into();
            self.ck.pp.commit(coeff_p_bar, w_bar).unwrap()
        };

        // compute challenge alpha, must be in unit group of field
        let mut alpha = Scalar::random(rng);
        while alpha == Scalar::zero() {
            alpha = Scalar::random(rng);
        }

        let coeff_p_bar: Vec<Scalar> = p_bar.clone().into();
        let alpha_coeff_p_poly: Vec<Scalar> = coeff_p_bar.iter().map(|s| alpha * s).collect();
        let p_prime_poly: Polynomial<Scalar> = p.clone() + Polynomial::from(alpha_coeff_p_poly);
        let w_prime = w + alpha * w_bar;
        // compute non hiding commitment to p:  C + α\hat{C} − ω'S
        let non_hiding_commitment_C_prime: RistrettoPoint = commitment.0 + (alpha * hiding_commitment_C_bar.0) - (w_prime * self.ck.pp.G.0);
        let zeroth_challenge = Scalar::random(rng);
        let H_prime = zeroth_challenge * self.ck.pp.H[0].0;
        // get coefficients of original polynomial p
        let coeffs_p_poly: Vec<Scalar> = p.into();
        // compute d+1 powers of z from 0 -> d
        let powers_of_z = {
            let mut temp: Vec<Scalar> = Vec::with_capacity((deg + 1) as usize);
            for i in 0..((deg + 1) as usize) {
                if i == 0 {
                    temp[i] = Scalar::zero();
                } else {
                    temp[i] = temp[i-1] * z;
                }
            }

            temp
        };
        let G_0 = self.ck.pp.H.clone();
        let log_2_d1 = log2(deg+1);
        // need to understand following steps from open construction
        Ok(())
    }
}

fn log2(num: u32) -> u8 {
    if num == 1 {
        0
    } else {
        1 + log2(num / 2)
    }
}
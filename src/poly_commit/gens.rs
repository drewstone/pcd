#![allow(non_snake_case)]

use super::errors::{ProofError, ProofResult};
use curve25519_dalek::constants;
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use codec::{Encode, Decode};
use sha3::Sha3_512;
use crate::types::GroupElt;
use polynomials::sparse::SparsePolynomial;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

/// A collection of generator points that can be used to compute various proofs
/// in this module. To create an instance of [`UniversalParams`] it is recommended to
/// call UniversalParams::new(`n`), where `n` is the number of bits to be used in
/// proofs and verifications.
#[derive(Debug, Clone, Encode, Decode)]
pub struct UniversalParams {
    pub n_bits: u32,
    pub D: u32,
    pub G: GroupElt,
    pub comm_key: Vec<GroupElt>,
    pub h: GroupElt,
    pub S: GroupElt,
}

impl UniversalParams {
    /// Create a new instance of [`UniversalParams`] with enough generator points to
    /// support proof and verification over an `n_bit` sized set.
    pub fn new(n_bits: u32, L: u32) -> ProofResult<UniversalParams> {
        if n_bits <= 1 {
            return Err(ProofError::SetIsTooSmall);
        }
        if n_bits > 32 {
            return Err(ProofError::SetIsTooLarge);
        };

        // Compute L sufficiently random generator points
        // r*G + msg[0]*H[0] + ... + msg[2n-1]*H[2n-1]
        //
        // G       = Ristretto Base Point
        // h    = hash(G)
        // H[0]    = hash(H[0])
        //  .           .
        //  .           .
        //  .           .
        // H[L-1] = hash(H[L-2])
        let h = GroupElt(RistrettoPoint::hash_from_bytes::<Sha3_512>(
            constants::RISTRETTO_BASEPOINT_POINT.compress().as_bytes(),
        ));
        let mut pts = Vec::new();
        let prev = h;
        for i in 0..(L as usize) {
            pts.push(GroupElt(RistrettoPoint::hash_from_bytes::<Sha3_512>(
                prev.0.compress().as_bytes(),
            )));
            prev = pts[pts.len() - 1];
        }

        let S = GroupElt(RistrettoPoint::hash_from_bytes::<Sha3_512>(
            prev.0.compress().as_bytes(),
        ));

        Ok(UniversalParams {
            n_bits,
            D: L,
            G: GroupElt(constants::RISTRETTO_BASEPOINT_POINT),
            h: h,
            comm_key: pts,
            S: S,
        })
    }

    /// Returns the maximum set size that can be processed in a proof or
    /// verification. For example, a 10 bit proof would only be able to support
    /// proofs over a set with at most `2^10 = 1024` members. Note, proofs over
    /// smaller sets will be extended by repeating the first member.
    pub fn max_set_size(&self) -> usize {
        2usize.checked_pow(self.n_bits as u32).unwrap()
    }

    pub fn setup(n_bits: u32, d: u32) -> ProofResult<UniversalParams> {
        Self::new(n_bits, d+1)
    }

    // (ckPC,rkPC) := (ck, H),(ck, H)
    pub fn trim(&mut self, l: u32) -> Self {
        // Always use pp_pc = (pp, H[0])
        assert!(l <= self.D);
        let mut temp = self.comm_key.clone();
        temp.truncate(l as usize);
        Self {
            n_bits: self.n_bits,
            D: self.D,
            G: self.G,
            h: self.h,
            S: self.S,
            comm_key: temp,
        }
    }

    pub fn get_full_coeffs(&self, p: SparsePolynomial<Scalar>) -> ProofResult<Vec<Scalar>> {
        let coeffs: Vec<Scalar> = Vec::with_capacity(p.degree());
        let elts: Vec<(u64, Scalar)> = p.into();
        for elt in elts.iter() {
            let (exp, c) = elt;
            coeffs[*exp as usize] = *c;
        }

        Ok(coeffs)
    }

    pub fn commit(&self, msg: Vec<Scalar>, r: Scalar) -> ProofResult<GroupElt> {
        assert!(msg.len() <= self.comm_key.len());
        // compute m_i * g_i for all messages
        let temp: Vec<GroupElt> = msg.iter()
            .enumerate()
            .map(|(inx, m_i)| GroupElt(m_i * self.comm_key[inx].0))
            .collect();
        // compute sum(m_i * g_i) for all i from before
        let summed: GroupElt = temp.iter().fold(GroupElt(Default::default()), |a, b| GroupElt(a.0 + b.0));
        // add blinding factor
        let hiding_sum: GroupElt = GroupElt(r * self.G.0 + summed.0);
        Ok(hiding_sum)
    }

    pub fn commit_with_points(&self, msg: Vec<Scalar>, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt> {
        assert!(msg.len() == pts.len());
        // compute m_i * g_i for all messages
        let temp: Vec<GroupElt> = msg.iter()
            .enumerate()
            .map(|(inx, m_i)| GroupElt(m_i * pts[inx].0))
            .collect();
        // compute sum(m_i * g_i) for all i from before
        let summed: GroupElt = temp.iter().fold(GroupElt(Default::default()), |a, b| GroupElt(a.0 + b.0));
        // add blinding factor
        let hiding_sum: GroupElt = GroupElt(r * self.G.0 + summed.0);
        Ok(hiding_sum)
    }

    pub fn commit_sparse(&self, p: SparsePolynomial<Scalar>, r: Scalar) -> ProofResult<GroupElt> {
        match self.get_full_coeffs(p) {
            Ok(coeffs) => self.commit(coeffs, r),
            _ => panic!("Unimplemented error handling"),
        }
    }

    pub fn commit_sparse_with_points(&self, p: SparsePolynomial<Scalar>, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt> {
        match self.get_full_coeffs(p) {
            Ok(coeffs) => self.commit_with_points(coeffs, r, pts),
            _ => panic!("Unimplemented error handling"),
        }
    }
}

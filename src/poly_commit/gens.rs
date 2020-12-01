#![allow(non_snake_case)]

use super::errors::{ProofError, ProofResult};
use curve25519_dalek::constants;
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use codec::{Encode, Decode};
use sha3::Sha3_512;
use crate::types::GroupElt;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

/// A collection of generator points that can be used to compute various proofs
/// in this module. To create an instance of [`ProofGens`] it is recommended to
/// call ProofGens::new(`n`), where `n` is the number of bits to be used in
/// proofs and verifications.
#[derive(Debug, Clone, Encode, Decode)]
pub struct ProofGens {
    pub n_bits: u32,
    pub D: u32,
    pub G: GroupElt,
    pub H: Vec<GroupElt>,
}

impl ProofGens {
    /// Create a new instance of [`ProofGens`] with enough generator points to
    /// support proof and verification over an `n_bit` sized set.
    ///
    /// ```
    /// # use one_of_many_proofs::proofs::ProofGens;
    /// // Support 10 bit membership proofs
    /// let gens = ProofGens::new(10);
    /// ```
    pub fn new(n_bits: u32, L: u32) -> ProofResult<ProofGens> {
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
        // H[0]    = hash(G)
        // H[1]    = hash(H[0])
        //  .           .
        //  .           .
        //  .           .
        // H[L] = hash(H[L-1])
        let mut pts = Vec::with_capacity(L as usize);
        pts.push(GroupElt(RistrettoPoint::hash_from_bytes::<Sha3_512>(
            constants::RISTRETTO_BASEPOINT_POINT.compress().as_bytes(),
        )));
        for i in 1..(L as usize) {
            pts.push(GroupElt(RistrettoPoint::hash_from_bytes::<Sha3_512>(
                pts[i - 1].0.compress().as_bytes(),
            )));
        }

        // pp = (Σ=(G_1,G_2,...,G_L)=H[1:], S=G, D=degree)
        Ok(ProofGens {
            n_bits,
            D: L,
            G: GroupElt(constants::RISTRETTO_BASEPOINT_POINT),
            H: pts,
        })
    }

    /// Returns the maximum set size that can be processed in a proof or
    /// verification. For example, a 10 bit proof would only be able to support
    /// proofs over a set with at most `2^10 = 1024` members. Note, proofs over
    /// smaller sets will be extended by repeating the first member.
    pub fn max_set_size(&self) -> usize {
        2usize.checked_pow(self.n_bits as u32).unwrap()
    }

    pub fn setup(n_bits: u32, d: u32) -> ProofResult<ProofGens> {
        Self::new(n_bits, d+1)
    }

    // (ckPC,rkPC) := (ck, H),(ck, H)
    pub fn trim(&mut self, l: u32) -> Self {
        // Always use pp_pc = (pp, H[0])
        assert!(l + 1 <= self.D);
        let mut temp = self.H.clone();
        temp.truncate((l + 1) as usize);
        Self {
            n_bits: self.n_bits,
            G: self.G,
            H: temp,
            D: self.D,
        }
    }

    pub fn commit(&self, msg: Vec<Scalar>, r: Scalar) -> ProofResult<GroupElt> {
        assert!(msg.len() <= self.H.len() - 1);
        // compute m_i * g_i for all messages
        let temp: Vec<GroupElt> = msg.iter()
            .enumerate()
            .map(|(inx, m_i)| GroupElt(m_i * self.H[inx + 1].0))
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
}

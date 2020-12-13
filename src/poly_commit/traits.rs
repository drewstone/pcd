use super::*;
use curve25519_dalek::ristretto::RistrettoPoint;
use crate::poly_commit::errors::{ProofResult};
use curve25519_dalek::scalar::Scalar;
use crate::types::GroupElt;
use polynomials::sparse::SparsePolynomial;

pub trait PolyCommitment {
    fn commit(&mut self, p: SparsePolynomial<Scalar>, deg: u32, r: Scalar) -> ProofResult<GroupElt>;
    fn commit_with_points(&mut self, p: SparsePolynomial<Scalar>, deg: u32, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt>;
    fn commit_vec_with_points(&mut self, v: Vec<Scalar>, deg: u32, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt>;
    fn challenge_scalar_i(&mut self, i: u8) -> Scalar;
    fn inner_product(&self, left: Vec<Scalar>, right: Vec<Scalar>) -> Scalar;
    fn check(&mut self, commitment: GroupElt, deg: u32, z: Scalar, v: Scalar, proof: EvaluationProof) -> bool;
    fn succinct_check(&mut self, commitment: GroupElt, deg: u32, z: Scalar, v: Scalar, proof: EvaluationProof) -> ProofResult<(SparsePolynomial<Scalar>, GroupElt)>;
    fn open<R: RngCore + CryptoRng>(
        &mut self,
        p_poly: SparsePolynomial<Scalar>,
        commitment: GroupElt,
        deg: u32,
        z: Scalar,
        w: Scalar,
        rng: &mut R
    ) -> ProofResult<EvaluationProof>;
}
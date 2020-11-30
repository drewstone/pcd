//! Defines a `TranscriptProtocol` trait for using a Merlin transcript.
use super::errors::{ProofError, ProofResult};

use curve25519_dalek::ristretto::CompressedRistretto;
use curve25519_dalek::scalar::Scalar;

use merlin::Transcript;

pub trait TranscriptProtocol {
    /// Append a domain separator for a poly_commit proof of degree d poly
    fn poly_commit_domain_sep(&mut self, d: u64);

    /// Append a `scalar` with the given `label`.
    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar);

    /// Append a `point` with the given `label`.
    fn append_point(&mut self, label: &'static [u8], point: &CompressedRistretto);

    /// Check that a point is not the identity, then append it to the
    /// transcript.  Otherwise, return an error.
    fn validate_and_append_point(
        &mut self,
        label: &'static [u8],
        point: &CompressedRistretto,
    ) -> ProofResult<()>;

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar;
}

impl TranscriptProtocol for Transcript {
    fn poly_commit_domain_sep(&mut self, d: u64) {
        self.append_message(b"dom-sep", b"poly_commit v1");
        self.append_u64(b"d", d);
    }

    fn append_scalar(&mut self, label: &'static [u8], scalar: &Scalar) {
        self.append_message(label, scalar.as_bytes());
    }

    fn append_point(&mut self, label: &'static [u8], point: &CompressedRistretto) {
        self.append_message(label, point.as_bytes());
    }

    fn validate_and_append_point(
        &mut self,
        label: &'static [u8],
        point: &CompressedRistretto,
    ) -> ProofResult<()> {
        use curve25519_dalek::traits::IsIdentity;

        if point.is_identity() {
            Err(ProofError::VerificationFailed)
        } else {
            Ok(self.append_message(label, point.as_bytes()))
        }
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 64];
        self.challenge_bytes(label, &mut buf);

        Scalar::from_bytes_mod_order_wide(&buf)
    }
}

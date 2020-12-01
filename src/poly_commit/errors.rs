use curve25519_dalek::scalar::Scalar;

pub type ProofResult<T> = core::result::Result<T, ProofError>;

#[derive(Debug, PartialEq, Clone)]
pub enum ProofError {
    InvalidProofSize,
    InvalidScalar(Scalar),
    IndexOutOfBounds,
    PointIsIdentity,
    SetIsTooLarge,
    SetIsTooSmall,
    Unimplemented,
    VerificationFailed,
}
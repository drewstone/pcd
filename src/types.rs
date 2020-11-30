use codec::{Decode, Encode, EncodeLike, Input};
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use curve25519_dalek::scalar::Scalar;

#[derive(Eq, PartialEq, Clone, Default, Debug, Copy)]
pub struct GroupElt(pub RistrettoPoint);
#[derive(Eq, PartialEq, Clone, Default, Debug, Copy)]
pub struct PrivateKey(pub Scalar);

pub const SIZE: usize = 32;

impl Encode for GroupElt {
    fn using_encoded<R, F: FnOnce(&[u8]) -> R>(&self, f: F) -> R {
        (self.0).compress().0.using_encoded(f)
    }
}

impl EncodeLike for GroupElt {}

impl Decode for GroupElt {
    fn decode<I: Input>(input: &mut I) -> Result<Self, codec::Error> {
        match <[u8; SIZE] as Decode>::decode(input).map(CompressedRistretto) {
            Ok(elt) => Ok(GroupElt(elt.decompress().unwrap())),
            Err(e) => Err(e),
        }
    }
}

impl Encode for PrivateKey {
    fn using_encoded<R, F: FnOnce(&[u8]) -> R>(&self, f: F) -> R {
        (self.0).as_bytes().using_encoded(f)
    }
}

impl EncodeLike for PrivateKey {}

impl Decode for PrivateKey {
    fn decode<I: Input>(input: &mut I) -> Result<Self, codec::Error> {
        match <[u8; SIZE] as Decode>::decode(input) {
            Ok(elt) => Ok(PrivateKey(
                Scalar::from_canonical_bytes(elt).unwrap_or(Scalar::zero()),
            )),
            Err(e) => Err(e),
        }
    }
}

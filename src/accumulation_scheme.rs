pub struct IVCProof<T: Trait> {
	acc: T::Accumulator,
	proof: T::Proof,
}

pub trait Trait {
	type Instance;
	type Proof;
	type Accumulator;
}

pub trait PCDFromAccumulationScheme: Trait {
	fn generate();
	fn prove(inst: Self::Instance, proof: Self::Proof) -> Self::Accumulator;
	fn verify(inst: Self::Instance, proof: Self::Proof, acc: Self::Accumulator, acc_star: Self::Accumulator) -> bool;
	fn decide(acc: Self::Accumulator) -> bool;

	fn pcd_prove(instance: Self::Instance, proof: Self::Proof, acc: Self::Accumulator) -> Self::Proof;
}

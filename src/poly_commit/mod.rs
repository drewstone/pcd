#![allow(non_snake_case)]

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use crate::types::GroupElt;
use polynomials::sparse::SparsePolynomial;
use merlin::Transcript;
use transcript::*;

pub mod traits;
pub mod errors;
pub mod transcript;
pub mod gens;

use traits::*;
use errors::*;

#[cfg(feature = "std")]
use rand::{CryptoRng, RngCore};

pub struct PolyCommitmentScheme {
    pp: gens::UniversalParams,
    transcript: Transcript,
}

pub struct EvaluationProof {
    L: Vec<GroupElt>,
    R: Vec<GroupElt>,
    U: GroupElt,
    c: Scalar,
    C: GroupElt,
    w_prime: Scalar, 
}

impl PolyCommitment for PolyCommitmentScheme {
    fn commit(&mut self, p: SparsePolynomial<Scalar>, deg: u32, r: Scalar) -> ProofResult<GroupElt> {
        assert!(p.degree() == deg as usize);
        let commit = self.pp.commit_sparse(p, r);
        self.transcript.append_point(b"poly_commit::commit_sparse_poly", &commit.clone().unwrap().0.compress());
        commit
    }

    fn commit_with_points(&mut self, p: SparsePolynomial<Scalar>, deg: u32, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt> {
        assert!(pts.len() <= deg as usize);
        assert!(p.degree() <= deg as usize);
        let commit = self.pp.commit_sparse_with_points(p, r, pts);
        self.transcript.append_point(b"poly_commit::commit_sparse_poly_with_points", &commit.clone().unwrap().0.compress());
        commit
    }

    fn commit_vec_with_points(&mut self, v: Vec<Scalar>, deg: u32, r: Scalar, pts: Vec<GroupElt>) -> ProofResult<GroupElt> {
        assert!(pts.len() <= deg as usize);
        assert!(v.len() <= deg as usize);
        assert!(v.len() == pts.len());
        let commit = self.pp.commit_with_points(v, r, pts);
        self.transcript.append_point(b"poly_commit::commit_direct_with_points", &commit.clone().unwrap().0.compress());
        commit
    }

    #[cfg(feature = "std")]
    fn open<R: RngCore + CryptoRng>(&mut self, p: SparsePolynomial<Scalar>, commitment: GroupElt, deg: u32, z: Scalar, w: Scalar, rng: &mut R) -> ProofResult<EvaluationProof> {
        // evaluation is only necessary for randomness sampling
        let _p_eval = p.eval(z);
        // sample random poly with p_bar(z) = 0
        let p_bar = {
            let mut partial_poly: SparsePolynomial<Scalar> = SparsePolynomial::new();
            for i in 0..((deg - 1) as usize) {
                partial_poly.push((i as u64, Scalar::random(rng)));
            }

            // z_root_polly = (x - z)
            let mut z_root_poly: SparsePolynomial<Scalar> = SparsePolynomial::new();
            z_root_poly.push((0u64, z));
            z_root_poly.push((1u64, Scalar::from(1u32)));

            partial_poly * z_root_poly
        };

        let w_bar = Scalar::random(rng);
        // compute hiding commitment: \hat{C} = CM.Commit(ck, coeff_p_bar, w_bar)
        let hiding_commitment_C_bar = self.commit(p_bar.clone(), deg, w_bar).unwrap();
        // compute challenge alpha, must be in unit group of field
        // α := ρ(C, z, v, C¯)
        let mut alpha = self.transcript.challenge_scalar(b"poly_commit::alpha");
        while alpha == Scalar::zero() {
            alpha = self.transcript.challenge_scalar(b"poly_commit::alpha");
        }

        let coeff_p_bar: Vec<(u64, Scalar)> = p_bar.clone().into();
        let alpha_coeff_p_poly: Vec<(u64, Scalar)> = coeff_p_bar.iter().map(|s| (s.0, alpha * s.1)).collect();
        let _p_prime_poly: SparsePolynomial<Scalar> = p.clone() + SparsePolynomial::from(alpha_coeff_p_poly);
        let w_prime = w + alpha * w_bar;
        // compute non hiding commitment to p:  C + α\hat{C} − ω'S
        let _non_hiding_commitment_C_prime: RistrettoPoint = commitment.0 + (alpha * hiding_commitment_C_bar.0) - (w_prime * self.pp.S.0);
        // ρ0(C', z, v) is the random oracle,
        // TODO: Use proper randomness as outlined in the paper
        let zeroth_challenge = self.challenge_scalar_i(0);
        let h_prime = zeroth_challenge * self.pp.h.0;
        // get coefficients of original SparsePolynomial p
        let mut coeffs_p_poly: Vec<Scalar> = self.pp.get_full_coeffs(p).unwrap();
        // compute d+1 powers of z from 0 -> d
        let mut powers_of_z = {
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

        let mut gs: Vec<GroupElt> = self.pp.comm_key.iter().map(|e| *e).collect();
        let mut left_commits: Vec<GroupElt> = vec![];
        let mut right_commits: Vec<GroupElt> = vec![];
        let log2_val = log2(deg + 1);
        // compute left and right commitments for proof
        for i in 0..log2_val {
            let left: Vec<GroupElt> = gs[1..gs.len() / 2]
                .iter()
                .chain(&[GroupElt(h_prime)])
                .map(|e| *e)
                .collect();
            let right: Vec<GroupElt> = gs[gs.len() / 2..]
                .iter()
                .map(|e| *e)
                .collect();

            let left_coeffs_p_poly: Vec<Scalar> = coeffs_p_poly[..coeffs_p_poly.len() / 2]
                .iter()
                .map(|e| *e)
                .collect();
            let right_coeffs_p_poly: Vec<Scalar> = coeffs_p_poly[coeffs_p_poly.len() / 2..]
                .iter()
                .map(|e| *e)
                .collect();

            let left_powers_of_z: Vec<Scalar> = powers_of_z[..powers_of_z.len() / 2]
                .iter()
                .map(|e| *e)
                .collect();
            let right_powers_of_z: Vec<Scalar> = powers_of_z[powers_of_z.len() / 2..]
                .iter()
                .map(|e| *e)
                .collect();



            let left_commit = self.commit_with_points(
                SparsePolynomial::from(
                    right_coeffs_p_poly.clone()
                        .into_iter()
                        .enumerate()
                        .filter(|(inx, c)| *c == Scalar::zero())
                        .map(|(inx, c)| (inx as u64, c))
                        .collect()
                ),
                deg,
                self.inner_product(right_coeffs_p_poly.clone(), left_powers_of_z),
                left.clone(),
            );

            let right_commit = self.commit_with_points(
                SparsePolynomial::from(
                    left_coeffs_p_poly.clone()
                        .into_iter()
                        .enumerate()
                        .filter(|(inx, c)| *c == Scalar::zero())
                        .map(|(inx, c)| (inx as u64, c))
                        .collect()
                ),
                deg,
                self.inner_product(left_coeffs_p_poly.clone(), right_powers_of_z),
                right.clone(),
            );

            left_commits.push(left_commit.unwrap());
            right_commits.push(right_commit.unwrap());
            // ξi:= ρ0(ξi−1, Li, R)
            let ith_challenge = self.challenge_scalar_i(i+1);

            gs = left.iter()
                .chain(right
                    .iter()
                    .map(|elt| GroupElt(ith_challenge * elt.0))
                    .collect::<Vec<GroupElt>>()
                    .iter())
                .map(|elt| *elt)
                .collect();
            coeffs_p_poly = left_coeffs_p_poly.iter()
                .chain(right_coeffs_p_poly
                    .iter()
                    .map(|elt| ith_challenge.invert() * elt)
                    .collect::<Vec<Scalar>>()
                    .iter())
                .map(|elt| *elt)
                .collect();
            powers_of_z = left_coeffs_p_poly.iter()
                .chain(right_coeffs_p_poly
                    .iter()
                    .map(|elt| ith_challenge * elt)
                    .collect::<Vec<Scalar>>()
                    .iter())
                .map(|elt| *elt)
                .collect();
        }
        
        let proof = EvaluationProof {
            L: left_commits,
            R: right_commits,
            U: gs[(log2_val - 1) as usize],
            c: coeffs_p_poly[(log2_val - 1) as usize],
            C: hiding_commitment_C_bar,
            w_prime: w_prime, 
        };
        
        Ok(proof)
    }

    fn check(&mut self, commitment: GroupElt, deg: u32, z: Scalar, v: Scalar, proof: EvaluationProof) -> bool {
        // H[0] considered H, H[1:] considered hk
        let d_prime = self.pp.comm_key.len() - 1;
        // // set rk := (<group>, S, H, d')
        // let rk: (GroupElt, GroupElt, GroupElt, usize) = (self.pp.G, self.pp.G, self.pp.h, d_prime);
        // check PC.SuccinctCheck(rk, C, d, z, v, \pi)
        assert!(d_prime == deg as usize);
        let (h, U) = self.succinct_check(commitment, deg, z, v, proof).unwrap();
        // check that U = CM.Commit(ck,coeffs_h_poly
        let computed_commitment = self.pp.commit_sparse(h, Scalar::zero());
        return computed_commitment.unwrap().0 == U.0;
    }

    fn succinct_check(&mut self, commitment: GroupElt, deg: u32, z: Scalar, v: Scalar, proof: EvaluationProof) -> ProofResult<(SparsePolynomial<Scalar>, GroupElt)> {
        let alpha = self.transcript.challenge_scalar(b"poly_commit::alpha");
        let non_hiding_commitment_c_prime = GroupElt(commitment.0 + (alpha * proof.C.0) - (proof.w_prime * self.pp.S.0));
        let zeroth_challenge = self.challenge_scalar_i(0);
        let H_ = GroupElt(zeroth_challenge * self.pp.h.0);
        let C_0 = GroupElt(non_hiding_commitment_c_prime.0 + (v * H_.0));
        let mut ith_commitment = C_0;
        let log2_val = log2(deg + 1);
        // compute left and right commitments for proof
        let mut challenges: Vec<Scalar> = vec![];
        for i in 0..log2_val {
            // ξi:= ρ0(ξi−1, Li, R)
            let ith_challenge = self.challenge_scalar_i(i+1);
            ith_commitment = GroupElt((ith_challenge.invert() * proof.L[i as usize].0) + ith_commitment.0);
            challenges.push(ith_challenge);
        }
        let mut rev_challenges = challenges.clone();
        rev_challenges.reverse();
        let mut last_power = 0;
        let h = rev_challenges.iter().enumerate().map(|(inx, ch)| {
            last_power = if last_power == 0 { 1 } else { 2 * last_power };
            SparsePolynomial::from(vec![(1,Scalar::zero()), (last_power, *ch)])
        }).fold(SparsePolynomial::from(vec![(1, Scalar::one())]), |total, next| { total * next });
        let v_prime = h.eval(z).unwrap();
        let last_commit = self.commit_vec_with_points(
            vec![proof.c, v_prime],
            2,
            Scalar::zero(),
            vec![proof.U, H_],
        ).unwrap();
        assert!(last_commit == ith_commitment);
        return Ok((h, proof.U))
    }

    fn inner_product(&self, left: Vec<Scalar>, right: Vec<Scalar>) -> Scalar {
        assert!(left.len() == right.len());

        let mut summed: Scalar = Scalar::zero();
        for i in 0..left.len() {
            summed += left[i] * right[i];
        }
        summed
    }

    fn challenge_scalar_i(&mut self, i: u8) -> Scalar {
        let s: &'static [u8] = match i {
            0 => b"poly_commit::challenge_0",
            1 => b"poly_commit::challenge_1",
            2 => b"poly_commit::challenge_2",
            3 => b"poly_commit::challenge_3",
            4 => b"poly_commit::challenge_4",
            5 => b"poly_commit::challenge_5",
            6 => b"poly_commit::challenge_6",
            7 => b"poly_commit::challenge_7",
            8 => b"poly_commit::challenge_8",
            9 => b"poly_commit::challenge_9",
            10 => b"poly_commit::challenge_10",
            11 => b"poly_commit::challenge_11",
            12 => b"poly_commit::challenge_12",
            13 => b"poly_commit::challenge_13",
            14 => b"poly_commit::challenge_14",
            15 => b"poly_commit::challenge_15",
            16 => b"poly_commit::challenge_16",
            17 => b"poly_commit::challenge_17",
            18 => b"poly_commit::challenge_18",
            19 => b"poly_commit::challenge_19",
            20 => b"poly_commit::challenge_20",
            21 => b"poly_commit::challenge_21",
            22 => b"poly_commit::challenge_22",
            23 => b"poly_commit::challenge_23",
            24 => b"poly_commit::challenge_24",
            25 => b"poly_commit::challenge_25",
            26 => b"poly_commit::challenge_26",
            27 => b"poly_commit::challenge_27",
            28 => b"poly_commit::challenge_28",
            29 => b"poly_commit::challenge_29",
            30 => b"poly_commit::challenge_30",
            31 => b"poly_commit::challenge_31",
            32 => b"poly_commit::challenge_32",
            33 => b"poly_commit::challenge_33",
            34 => b"poly_commit::challenge_34",
            35 => b"poly_commit::challenge_35",
            36 => b"poly_commit::challenge_36",
            37 => b"poly_commit::challenge_37",
            38 => b"poly_commit::challenge_38",
            39 => b"poly_commit::challenge_39",
            40 => b"poly_commit::challenge_40",
            41 => b"poly_commit::challenge_41",
            42 => b"poly_commit::challenge_42",
            43 => b"poly_commit::challenge_43",
            44 => b"poly_commit::challenge_44",
            45 => b"poly_commit::challenge_45",
            46 => b"poly_commit::challenge_46",
            47 => b"poly_commit::challenge_47",
            48 => b"poly_commit::challenge_48",
            49 => b"poly_commit::challenge_49",
            50 => b"poly_commit::challenge_50",
            51 => b"poly_commit::challenge_51",
            52 => b"poly_commit::challenge_52",
            53 => b"poly_commit::challenge_53",
            54 => b"poly_commit::challenge_54",
            55 => b"poly_commit::challenge_55",
            56 => b"poly_commit::challenge_56",
            57 => b"poly_commit::challenge_57",
            58 => b"poly_commit::challenge_58",
            59 => b"poly_commit::challenge_59",
            60 => b"poly_commit::challenge_60",
            61 => b"poly_commit::challenge_61",
            62 => b"poly_commit::challenge_62",
            63 => b"poly_commit::challenge_63",
            _ => panic!(""),
        };

        self.transcript.challenge_scalar(s)
    }
}

fn log2(num: u32) -> u8 {
    if num == 1 {
        0
    } else {
        1 + log2(num / 2)
    }
}
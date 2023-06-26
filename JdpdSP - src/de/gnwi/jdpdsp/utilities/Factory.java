/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2023  Achim Zielesny (achim.zielesny@googlemail.com)
 * 
 * Source code is available at <https://github.com/zielesny/Jdpd>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.gnwi.jdpdsp.utilities;


import de.gnwi.jdpdsp.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpdsp.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IOutput;
import de.gnwi.jdpdsp.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionPnhlnCalculator;
import de.gnwi.jdpdsp.interfaces.IRandom;
import de.gnwi.jdpdsp.interfaces.ITimeStepCalculator;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.particlePosition.ParticlePositionPool;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondForceConservativeCalculator;
import de.gnwi.jdpdsp.samples.integrationType.GwmvvTimeStepCalculator;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondPotentialCalculator;
import de.gnwi.jdpdsp.samples.integrationType.PnhlnTimeStepCalculator;
import de.gnwi.jdpdsp.samples.integrationType.S1mvvTimeStepCalculator;
import de.gnwi.jdpdsp.samples.integrationType.ScmvvTimeStepCalculator;
import de.gnwi.jdpdsp.samples.interactions.electrostatics.ParticlePairElectrostaticsAdHocForceConservativeCalculator;
import de.gnwi.jdpdsp.samples.interactions.electrostatics.ParticlePairElectrostaticsAdHocPotentialCalculator;
import de.gnwi.jdpdsp.samples.interactions.electrostatics.ParticlePairElectrostaticsDpdForceConservativeCalculator;
import de.gnwi.jdpdsp.samples.interactions.electrostatics.ParticlePairElectrostaticsDpdPotentialCalculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairDpdForceConservativeCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairDpdForceDissipativeCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairDpdForceRandomCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairScmvvDpdForceDissipativeCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairGwmvvDpdForceFullCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairScmvvDpdForceFullCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairDpdPotentialCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.dpdCutoff1.ParticlePairS1mvvVelocityUpdateCutoff1Calculator;
import de.gnwi.jdpdsp.samples.interactions.nearestNeighbor.ParticlePairNearestNeighborCalculator;
import de.gnwi.jdpdsp.samples.random.ApacheCommonsRandom;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.rng.simple.RandomSource;

/**
 * Factory for specific object instance creation
 * 
 * @author Achim Zielesny
 */
public class Factory {
    
    // <editor-fold defaultstate="collapsed" desc="Public Enums">
    /**
     * Random number generator type
     */
    public enum RandomType {

        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_ISAAC,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_JDK,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_JSF_32,            // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_JSF_64,            // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_KISS,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_MT,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_MT_64,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_MSWS,              // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_MWC_256,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_PCG_XSH_RR_32,     // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_PCG_XSH_RS_32,     // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_PCG_RXS_M_XS_64,   // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_PCG_MCG_XSH_RR_32, // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_PCG_MCG_XSH_RS_32, // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_SFC_32,            // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_SFC_64,            // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_SPLIT_MIX_64,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_TWO_CMRES,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_1024_A,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_19937_A,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_19937_C,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_44497_A,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_44497_B,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_WELL_512_A,
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XOR_SHIFT_1024_S_PHI,  // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_64_S,     // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_64_SS,    // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_128_PLUS,    // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_128_SS,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_128_PLUS, // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_128_SS,   // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_256_PLUS,    // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_256_SS,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_512_PLUS,    // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_512_SS,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_128_PP,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_128_PP,   // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_256_PP,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_SHI_RO_512_PP,      // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_1024_PP,  // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_1024_S,   // new in ACRNG 1.3
        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_XO_RO_SHI_RO_1024_SS;  // new in ACRNG 1.3

        /**
         * Random number generator type representations
         * 
         * @return Random number generator type representations
         */
        public static String[] getRandomNumberGeneratorTypeRepresentations() {
            String[] tmpRandomTypeRepresentations = new String[RandomType.values().length];
            int tmpIndex = 0;
            for (RandomType tmpRandomType : RandomType.values()) {
                tmpRandomTypeRepresentations[tmpIndex++] = tmpRandomType.toString();
            }
            Arrays.sort(tmpRandomTypeRepresentations);
            return tmpRandomTypeRepresentations;
        }
        
        /**
         * Returns string array of length 3 with non-jumpable RNGs (index 0), 
         * jumpable RNGs (index 1) and long-jumpable RNGs (index 2)
         * 
         * @return String array of length 3 with non-jumpable RNGs (index 0), 
         * jumpable RNGs (index 1) and long-jumpable RNGs (index 2)
         */
        public static String[] getJumpableRandomNumberGeneratorInfo() {
            StringBuilder tmpNonJumpableBuffer = new StringBuilder(1000);
            StringBuilder tmpJumpableBuffer = new StringBuilder(1000);
            StringBuilder tmpLongJumpableBuffer = new StringBuilder(1000);
            for (RandomType tmpRandomType : RandomType.values()) {
                boolean tmpIsJumpable = false;
                boolean tmpIsLongJumpable = false;
                switch(tmpRandomType) {
                    case ACRNG_ISAAC:
                        tmpIsJumpable = RandomSource.ISAAC.isJumpable();
                        tmpIsLongJumpable = RandomSource.ISAAC.isLongJumpable();
                        break;
                    case ACRNG_JDK:
                        tmpIsJumpable = RandomSource.JDK.isJumpable();
                        tmpIsLongJumpable = RandomSource.JDK.isLongJumpable();
                        break;
                    case ACRNG_JSF_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.JSF_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.JSF_32.isLongJumpable();
                        break;
                    case ACRNG_JSF_64: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.JSF_64.isJumpable();
                        tmpIsLongJumpable = RandomSource.JSF_64.isLongJumpable();
                        break;
                    case ACRNG_KISS:
                        tmpIsJumpable = RandomSource.KISS.isJumpable();
                        tmpIsLongJumpable = RandomSource.KISS.isLongJumpable();
                        break;
                    case ACRNG_MT:
                        tmpIsJumpable = RandomSource.MT.isJumpable();
                        tmpIsLongJumpable = RandomSource.MT.isLongJumpable();
                        break;
                    case ACRNG_MT_64:
                        tmpIsJumpable = RandomSource.MT_64.isJumpable();
                        tmpIsLongJumpable = RandomSource.MT_64.isLongJumpable();
                        break;
                    case ACRNG_MSWS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.MSWS.isJumpable();
                        tmpIsLongJumpable = RandomSource.MSWS.isLongJumpable();
                        break;
                    case ACRNG_MWC_256:
                        tmpIsJumpable = RandomSource.MWC_256.isJumpable();
                        tmpIsLongJumpable = RandomSource.MWC_256.isLongJumpable();
                        break;
                    case ACRNG_PCG_XSH_RR_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.PCG_XSH_RR_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.PCG_XSH_RR_32.isLongJumpable();
                        break;
                    case ACRNG_PCG_XSH_RS_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.PCG_XSH_RS_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.PCG_XSH_RS_32.isLongJumpable();
                        break;
                    case ACRNG_PCG_RXS_M_XS_64: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.PCG_RXS_M_XS_64.isJumpable();
                        tmpIsLongJumpable = RandomSource.PCG_RXS_M_XS_64.isLongJumpable();
                        break;
                    case ACRNG_PCG_MCG_XSH_RR_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.PCG_MCG_XSH_RR_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.PCG_MCG_XSH_RR_32.isLongJumpable();
                        break;
                    case ACRNG_PCG_MCG_XSH_RS_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.PCG_MCG_XSH_RS_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.PCG_MCG_XSH_RS_32.isLongJumpable();
                        break;
                    case ACRNG_SFC_32: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.SFC_32.isJumpable();
                        tmpIsLongJumpable = RandomSource.SFC_32.isLongJumpable();
                        break;
                    case ACRNG_SFC_64: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.SFC_64.isJumpable();
                        tmpIsLongJumpable = RandomSource.SFC_64.isLongJumpable();
                        break;
                    case ACRNG_SPLIT_MIX_64:
                        tmpIsJumpable = RandomSource.SPLIT_MIX_64.isJumpable();
                        tmpIsLongJumpable = RandomSource.SPLIT_MIX_64.isLongJumpable();
                        break;
                    case ACRNG_TWO_CMRES:
                        tmpIsJumpable = RandomSource.TWO_CMRES.isJumpable();
                        tmpIsLongJumpable = RandomSource.TWO_CMRES.isLongJumpable();
                        break;
                    case ACRNG_WELL_1024_A:
                        tmpIsJumpable = RandomSource.WELL_1024_A.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_1024_A.isLongJumpable();
                        break;
                    case ACRNG_WELL_19937_A:
                        tmpIsJumpable = RandomSource.WELL_19937_A.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_19937_A.isLongJumpable();
                        break;
                    case ACRNG_WELL_19937_C:
                        tmpIsJumpable = RandomSource.WELL_19937_C.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_19937_C.isLongJumpable();
                        break;
                    case ACRNG_WELL_44497_A:
                        tmpIsJumpable = RandomSource.WELL_44497_A.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_44497_A.isLongJumpable();
                        break;
                    case ACRNG_WELL_44497_B:
                        tmpIsJumpable = RandomSource.WELL_44497_B.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_44497_B.isLongJumpable();
                        break;
                    case ACRNG_WELL_512_A:
                        tmpIsJumpable = RandomSource.WELL_512_A.isJumpable();
                        tmpIsLongJumpable = RandomSource.WELL_512_A.isLongJumpable();
                        break;
                    case ACRNG_XOR_SHIFT_1024_S_PHI: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XOR_SHIFT_1024_S_PHI.isJumpable();
                        tmpIsLongJumpable = RandomSource.XOR_SHIFT_1024_S_PHI.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_64_S: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_64_S.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_64_S.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_64_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_64_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_64_SS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_128_PLUS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_128_PLUS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_128_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_128_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_128_SS.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_128_PLUS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_128_PLUS.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_128_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_128_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_128_SS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_256_PLUS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_256_PLUS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_256_PLUS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_256_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_256_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_256_SS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_512_PLUS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_512_PLUS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_512_PLUS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_512_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_512_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_512_SS.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_128_PP: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_128_PP.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_128_PP.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_128_PP: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_128_PP.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_128_PP.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_256_PP: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_256_PP.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_256_PP.isLongJumpable();
                        break;
                    case ACRNG_XO_SHI_RO_512_PP: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_SHI_RO_512_PP.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_SHI_RO_512_PP.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_1024_PP: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_1024_PP.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_1024_PP.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_1024_S: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_1024_S.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_1024_S.isLongJumpable();
                        break;
                    case ACRNG_XO_RO_SHI_RO_1024_SS: // new in ACRNG 1.3
                        tmpIsJumpable = RandomSource.XO_RO_SHI_RO_1024_SS.isJumpable();
                        tmpIsLongJumpable = RandomSource.XO_RO_SHI_RO_1024_SS.isLongJumpable();
                        break;
                    default:
                        throw new IllegalArgumentException("RandomType.getJumpableRandomNumberGenerators: Unknown random type.");
                }
                if (tmpIsLongJumpable) {
                    if (tmpLongJumpableBuffer.length() > 0) {
                        tmpLongJumpableBuffer.append(", ");
                    };
                    tmpLongJumpableBuffer.append(tmpRandomType.toString());
                } else if (tmpIsJumpable) {
                    if (tmpJumpableBuffer.length() > 0) {
                        tmpJumpableBuffer.append(", ");
                    };
                    tmpJumpableBuffer.append(tmpRandomType.toString());
                } else {
                    if (tmpNonJumpableBuffer.length() > 0) {
                        tmpNonJumpableBuffer.append(", ");
                    };
                    tmpNonJumpableBuffer.append(tmpRandomType.toString());
                }
            }
            return new String[] {
                tmpNonJumpableBuffer.toString(),
                tmpJumpableBuffer.toString(),
                tmpLongJumpableBuffer.toString()
            };
        }

        /**
         * Returns if random type is jumpable
         * 
         * @param aRandomType Random type
         * @return True: Random type is jumpable, false: Otherwise
         */
        public static boolean isJumpable(RandomType aRandomType) {
            switch(aRandomType) {
                case ACRNG_ISAAC:
                    return RandomSource.ISAAC.isJumpable() || RandomSource.ISAAC.isLongJumpable();
                case ACRNG_JDK:
                    return RandomSource.JDK.isJumpable() || RandomSource.JDK.isLongJumpable();
                case ACRNG_JSF_32: // new in ACRNG 1.3
                    return RandomSource.JSF_32.isJumpable() || RandomSource.JSF_32.isLongJumpable();
                case ACRNG_JSF_64: // new in ACRNG 1.3
                    return RandomSource.JSF_64.isJumpable() || RandomSource.JSF_64.isLongJumpable();
                case ACRNG_KISS:
                    return RandomSource.KISS.isJumpable() || RandomSource.KISS.isLongJumpable();
                case ACRNG_MT:
                    return RandomSource.MT.isJumpable() || RandomSource.MT.isLongJumpable();
                case ACRNG_MT_64:
                    return RandomSource.MT_64.isJumpable() || RandomSource.MT_64.isLongJumpable();
                case ACRNG_MSWS: // new in ACRNG 1.3
                    return RandomSource.MSWS.isJumpable() || RandomSource.MSWS.isLongJumpable();
                case ACRNG_MWC_256:
                    return RandomSource.MWC_256.isJumpable() || RandomSource.MWC_256.isLongJumpable();
                case ACRNG_PCG_XSH_RR_32: // new in ACRNG 1.3
                    return RandomSource.PCG_XSH_RR_32.isJumpable() || RandomSource.PCG_XSH_RR_32.isLongJumpable();
                case ACRNG_PCG_XSH_RS_32: // new in ACRNG 1.3
                    return RandomSource.PCG_XSH_RS_32.isJumpable() || RandomSource.PCG_XSH_RS_32.isLongJumpable();
                case ACRNG_PCG_RXS_M_XS_64: // new in ACRNG 1.3
                    return RandomSource.PCG_RXS_M_XS_64.isJumpable() || RandomSource.PCG_RXS_M_XS_64.isLongJumpable();
                case ACRNG_PCG_MCG_XSH_RR_32: // new in ACRNG 1.3
                    return RandomSource.PCG_MCG_XSH_RR_32.isJumpable() || RandomSource.PCG_MCG_XSH_RR_32.isLongJumpable();
                case ACRNG_PCG_MCG_XSH_RS_32: // new in ACRNG 1.3
                    return RandomSource.PCG_MCG_XSH_RS_32.isJumpable() || RandomSource.PCG_MCG_XSH_RS_32.isLongJumpable();
                case ACRNG_SFC_32: // new in ACRNG 1.3
                    return RandomSource.SFC_32.isJumpable() || RandomSource.SFC_32.isLongJumpable();
                case ACRNG_SFC_64: // new in ACRNG 1.3
                    return RandomSource.SFC_64.isJumpable() || RandomSource.SFC_64.isLongJumpable();
                case ACRNG_SPLIT_MIX_64:
                    return RandomSource.SPLIT_MIX_64.isJumpable() || RandomSource.SPLIT_MIX_64.isLongJumpable();
                case ACRNG_TWO_CMRES:
                    return RandomSource.TWO_CMRES.isJumpable() || RandomSource.TWO_CMRES.isLongJumpable();
                case ACRNG_WELL_1024_A:
                    return RandomSource.WELL_1024_A.isJumpable() || RandomSource.WELL_1024_A.isLongJumpable();
                case ACRNG_WELL_19937_A:
                    return RandomSource.WELL_19937_A.isJumpable() || RandomSource.WELL_19937_A.isLongJumpable();
                case ACRNG_WELL_19937_C:
                    return RandomSource.WELL_19937_C.isJumpable() || RandomSource.WELL_19937_C.isLongJumpable();
                case ACRNG_WELL_44497_A:
                    return RandomSource.WELL_44497_A.isJumpable() || RandomSource.WELL_44497_A.isLongJumpable();
                case ACRNG_WELL_44497_B:
                    return RandomSource.WELL_44497_B.isJumpable() || RandomSource.WELL_44497_B.isLongJumpable();
                case ACRNG_WELL_512_A:
                    return RandomSource.WELL_512_A.isJumpable() || RandomSource.WELL_512_A.isLongJumpable();
                case ACRNG_XOR_SHIFT_1024_S_PHI: // new in ACRNG 1.3
                    return RandomSource.XOR_SHIFT_1024_S_PHI.isJumpable() || RandomSource.XOR_SHIFT_1024_S_PHI.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_64_S: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_64_S.isJumpable() || RandomSource.XO_RO_SHI_RO_64_S.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_64_SS: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_64_SS.isJumpable() || RandomSource.XO_RO_SHI_RO_64_SS.isLongJumpable();
                case ACRNG_XO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_128_PLUS.isJumpable() || RandomSource.XO_SHI_RO_128_PLUS.isLongJumpable();
                case ACRNG_XO_SHI_RO_128_SS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_128_SS.isJumpable() || RandomSource.XO_SHI_RO_128_SS.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_128_PLUS.isJumpable() || RandomSource.XO_RO_SHI_RO_128_PLUS.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_128_SS: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_128_SS.isJumpable() || RandomSource.XO_RO_SHI_RO_128_SS.isLongJumpable();
                case ACRNG_XO_SHI_RO_256_PLUS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_256_PLUS.isJumpable() || RandomSource.XO_SHI_RO_256_PLUS.isLongJumpable();
                case ACRNG_XO_SHI_RO_256_SS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_256_SS.isJumpable() || RandomSource.XO_SHI_RO_256_SS.isLongJumpable();
                case ACRNG_XO_SHI_RO_512_PLUS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_512_PLUS.isJumpable() || RandomSource.XO_SHI_RO_512_PLUS.isLongJumpable();
                case ACRNG_XO_SHI_RO_512_SS: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_512_SS.isJumpable() || RandomSource.XO_SHI_RO_512_SS.isLongJumpable();
                case ACRNG_XO_SHI_RO_128_PP: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_128_PP.isJumpable() || RandomSource.XO_SHI_RO_128_PP.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_128_PP: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_128_PP.isJumpable() || RandomSource.XO_RO_SHI_RO_128_PP.isLongJumpable();
                case ACRNG_XO_SHI_RO_256_PP: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_256_PP.isJumpable() || RandomSource.XO_SHI_RO_256_PP.isLongJumpable();
                case ACRNG_XO_SHI_RO_512_PP: // new in ACRNG 1.3
                    return RandomSource.XO_SHI_RO_512_PP.isJumpable() || RandomSource.XO_SHI_RO_512_PP.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_1024_PP: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_1024_PP.isJumpable() || RandomSource.XO_RO_SHI_RO_1024_PP.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_1024_S: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_1024_S.isJumpable() || RandomSource.XO_RO_SHI_RO_1024_S.isLongJumpable();
                case ACRNG_XO_RO_SHI_RO_1024_SS: // new in ACRNG 1.3
                    return RandomSource.XO_RO_SHI_RO_1024_SS.isJumpable() || RandomSource.XO_RO_SHI_RO_1024_SS.isLongJumpable();
                default:
                    throw new IllegalArgumentException("RandomType.getJumpableRandomNumberGenerators: Unknown random type.");
            }
        }
        
        /**
         * Default random number generator type representation
         * 
         * @return Default random number generator type representation
         */
        public static String getDefaultRandomNumberGeneratorTypeRepresentation() {
            return Factory.RandomType.getDefaultRandomNumberGeneratorType().toString();
        }

        /**
         * Default random number generator type
         * 
         * @return Default random number generator type
         */
        public static RandomType getDefaultRandomNumberGeneratorType() {
            return RandomType.ACRNG_XO_SHI_RO_256_SS;
        }

        /**
         * Checks if aRandomNumberGeneratorTypeRepresentation is defined
         * 
         * @param aRandomNumberGeneratorTypeRepresentation Random number generator type representation
         * @return True: aRandomNumberGeneratorTypeRepresentation is defined, 
         * false: Otherwise
         */
        public static boolean isDefinedRandomNumberGeneratorTypeRepresentation(String aRandomNumberGeneratorTypeRepresentation) {
            // <editor-fold defaultstate="collapsed" desc="Checks">
            if (aRandomNumberGeneratorTypeRepresentation == null || aRandomNumberGeneratorTypeRepresentation.isEmpty()) {
                return false;
            }
            // </editor-fold>
            for (RandomType tmpRandomType : RandomType.values()) {
                if (aRandomNumberGeneratorTypeRepresentation.equals(tmpRandomType.toString())) {
                    return true;
                }
            }
            return false;
        }

    }  

    /**
     * DPD type
     */
    public enum DpdType {
        
        /**
         * DPD with cut-off length of 1.0
         */
        CUTOFF_LENGTH_ONE
        
    }
    
    /**
     * Electrostatics type
     */
    public enum ElectrostaticsType {
        
        /**
         * Ad-hoc electrostatics
         */
        AD_HOC,
        /**
         * DPD electrostatics
         */
        DPD
        
    }
    
    /**
     * Bond type
     */
    public enum BondType {
        
        /**
         * Harmonic bond
         */
        HARMONIC
        
    }

    /**
     * Integration type
     */
    public enum IntegrationType {
            
        /**
         * Groot-Warren Modified Velocity-Verlet (GWMVV)
         */
        GWMVV,
        /**
         * Self-consistent Modified Velocity-Verlet (SCMVV)
         */
        SCMVV,
        /**
         * Shardlow S1 Modified Velocity-Verlet (S1MVV)
         */
        S1MVV,
        /**
         * Nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN)
         */
        PNHLN;

        /**
         * Integration type representations
         * 
         * @return Integration type representations
         */
        public static String[] getIntegrationTypeRepresentations() {
            String[] tmpIntegrationTypeRepresentations = new String[IntegrationType.values().length];
            int tmpIndex = 0;
            for (IntegrationType tmIntegrationType : IntegrationType.values()) {
                tmpIntegrationTypeRepresentations[tmpIndex++] = tmIntegrationType.toString();
            }
            Arrays.sort(tmpIntegrationTypeRepresentations);
            return tmpIntegrationTypeRepresentations;
        }

        
        /**
         * Default integration type representation
         * 
         * @return Default integration type representation
         */
        public static String getDefaultIntegrationTypeRepresentation() {
            return Factory.IntegrationType.GWMVV.toString();
        }
        
    }    
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Random number generator type
     */
    private final RandomType randomType;

    /**
     * Number of random number generator warm-up steps
     */
    private final int numberORandomNumberGeneratorfWarmUpSteps;
    
    /**
     * DPD type
     */
    private final DpdType dpdType;

    /**
     * Electrostatics type
     */
    private final ElectrostaticsType electrostaticsType;
    
    /**
     * Bond type
     */
    private final BondType bondType;
    
    /**
     * Integration type
     */
    private final IntegrationType integrationType;
    
    /**
     * Integration parameters
     */
    private final Object[] integrationParameters;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Jumpable RNG
     */
    private ApacheCommonsRandom jumpableRandomNumberGenerator;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aRandomType Random number generator type
     * @param aNumberORandomNumberGeneratorfWarmUpSteps Number of random number 
     * generator warm-up steps (greater/equal 0)
     * @param aDpdType DPD type
     * @param anElectrostaticsType Electrostatics type
     * @param aBondType Bond type
     * @param anIntegrationType Integration type
     * @param anIntegrationParameters Integration parameters
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public Factory(
        RandomType aRandomType,
        int aNumberORandomNumberGeneratorfWarmUpSteps,
        DpdType aDpdType,
        ElectrostaticsType anElectrostaticsType,
        BondType aBondType,
        IntegrationType anIntegrationType,
        Object[] anIntegrationParameters
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aRandomType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aRandomType is null.");
        }
        if (aNumberORandomNumberGeneratorfWarmUpSteps < 0) {
            throw new IllegalArgumentException("Factory.Constructor: aNumberORandomNumberGeneratorfWarmUpSteps < 0.");
        }
        if (aDpdType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aDpdType is null.");
        }
        if (anElectrostaticsType == null) {
            throw new IllegalArgumentException("Factory.Constructor: anElectrostaticsType is null.");
        }
        if (aBondType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aBondType is null.");
        }
        if (anIntegrationType == null) {
            throw new IllegalArgumentException("Factory.Constructor: anIntegrationType is null.");
        }
        if (anIntegrationType == IntegrationType.GWMVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.GWMVV.");
            }
            if (anIntegrationParameters.length < 1) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 1 for IntegrationType.GWMVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Float)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Float for IntegrationType.GWMVV.");
            }
        }
        if (anIntegrationType == IntegrationType.SCMVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.SCMVV.");
            }
            if (anIntegrationParameters.length < 2) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 2 for IntegrationType.SCMVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Integer)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Integer for IntegrationType.SCMVV.");
            }
            if (!(anIntegrationParameters[1] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[1] is NOT instance of Boolean for IntegrationType.SCMVV.");
            }
        }
        if (anIntegrationType == IntegrationType.S1MVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.S1MVV.");
            }
            if (anIntegrationParameters.length < 1) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 1 for IntegrationType.S1MVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Boolean for IntegrationType.S1MVV.");
            }
        }
        if (anIntegrationType == IntegrationType.PNHLN) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.PNHLN.");
            }
            if (anIntegrationParameters.length < 2) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 2 for IntegrationType.PNHLN.");
            }
            if (!(anIntegrationParameters[0] instanceof Float)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Float for IntegrationType.PNHLN.");
            }
            if (!(anIntegrationParameters[1] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[1] is NOT instance of Boolean for IntegrationType.PNHLN.");
            }
        }
        // </editor-fold>
        this.randomType = aRandomType;
        this.numberORandomNumberGeneratorfWarmUpSteps = aNumberORandomNumberGeneratorfWarmUpSteps;
        this.dpdType = aDpdType;
        this.electrostaticsType = anElectrostaticsType;
        this.bondType = aBondType;
        this.integrationType = anIntegrationType;
        this.integrationParameters = anIntegrationParameters;
        this.jumpableRandomNumberGenerator = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    // <editor-fold defaultstate="collapsed" desc="- Random number generators">
    /**
     * Returns new or jumped random number generator instance
     * 
     * @param aSeed Seed value (greater/equal 0, is ignored for jumped random 
     * number generator)
     * @return New or jumped random number generator instance
     */
    public IRandom getNewOrJumpedRandomNumberGenerator(int aSeed) {
        if (this.jumpableRandomNumberGenerator != null) {
            return this.jumpableRandomNumberGenerator.getJumpedRng();
        } else {
            ApacheCommonsRandom tmpRng;
            switch(this.randomType) {
                case ACRNG_ISAAC:
                    tmpRng = new ApacheCommonsRandom(RandomSource.ISAAC, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_JDK:
                    tmpRng = new ApacheCommonsRandom(RandomSource.JDK, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_JSF_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.JSF_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_JSF_64: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.JSF_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_KISS:
                    tmpRng = new ApacheCommonsRandom(RandomSource.KISS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_MT:
                    tmpRng = new ApacheCommonsRandom(RandomSource.MT, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_MT_64:
                    tmpRng = new ApacheCommonsRandom(RandomSource.MT_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_MSWS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.MSWS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_MWC_256:
                    tmpRng = new ApacheCommonsRandom(RandomSource.MWC_256, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_PCG_XSH_RR_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.PCG_XSH_RR_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_PCG_XSH_RS_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.PCG_XSH_RS_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_PCG_RXS_M_XS_64: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.PCG_RXS_M_XS_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_PCG_MCG_XSH_RR_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.PCG_MCG_XSH_RR_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_PCG_MCG_XSH_RS_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.PCG_MCG_XSH_RS_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_SFC_32: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.SFC_32, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_SFC_64: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.SFC_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_SPLIT_MIX_64:
                    tmpRng = new ApacheCommonsRandom(RandomSource.SPLIT_MIX_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_TWO_CMRES:
                    tmpRng = new ApacheCommonsRandom(RandomSource.TWO_CMRES, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_1024_A:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_1024_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_19937_A:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_19937_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_19937_C:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_19937_C, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_44497_A:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_44497_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_44497_B:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_44497_B, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_WELL_512_A:
                    tmpRng = new ApacheCommonsRandom(RandomSource.WELL_512_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XOR_SHIFT_1024_S_PHI: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XOR_SHIFT_1024_S_PHI, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_64_S: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_64_S, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_64_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_64_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_128_PLUS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_128_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_128_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_128_PLUS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_128_PLUS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_128_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_128_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_256_PLUS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_256_PLUS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_256_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_256_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_512_PLUS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_512_PLUS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_512_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_512_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_128_PP: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_128_PP, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_128_PP: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_128_PP, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_256_PP: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_256_PP, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_SHI_RO_512_PP: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_SHI_RO_512_PP, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_1024_PP: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_1024_PP, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_1024_S: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_1024_S, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                case ACRNG_XO_RO_SHI_RO_1024_SS: // new in ACRNG 1.3
                    tmpRng = new ApacheCommonsRandom(RandomSource.XO_RO_SHI_RO_1024_SS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
                    break;
                default:
                    throw new IllegalArgumentException("Factory.getNewRandomNumberGenerator: Unknown random type.");
            }
            if (RandomType.isJumpable(this.randomType)) {
                this.jumpableRandomNumberGenerator = tmpRng;
            }
            return tmpRng;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- DPD calculators">
    // <editor-fold defaultstate="collapsed" desc="-- PNHLN calculators">
    /**
     * Returns particle pair velocity update plus G calculator for PNHLN integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair velocity update plus G calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionPnhlnCalculator getParticlePairPnhlnVelocityUpdatePlusGCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairPnhlnVelocityUpdatePlusGCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair velocity update plus G calculator for PNHLN integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair velocity update plus G calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionPnhlnCalculator getParticlePairPnhlnVelocityUpdatePlusGCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairPnhlnVelocityUpdatePlusGCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- S1MVV calculators">
    /**
     * Returns particle pair velocity update calculator for S1MVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair velocity update calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairS1mvvVelocityUpdateCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairS1mvvVelocityUpdateCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairS1mvvVelocityUpdateCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair velocity update calculator for S1MVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair velocity update calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairS1mvvVelocityUpdateCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairS1mvvVelocityUpdateCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairS1mvvVelocityUpdateCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- SCMVV calculators">
    /**
     * Returns particle pair DPD full force calculator for SCMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceFullCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceFullCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD full force calculator for SCMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceFullCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceFullCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD dissipative force calculator for 
     * SCMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceDissipativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD dissipative force calculator for 
     * SCMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceDissipativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- GWMVV calculators">
    /**
     * Returns particle pair DPD full force calculator for GWMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairGwmvvDpdForceFullCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairGwmvvDpdForceFullCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairGwmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD full force calculator for GWMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairGwmvvDpdForceFullCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairGwmvvDpdForceFullCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairGwmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- General calculators">
    /**
     * Returns particle pair DPD conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceConservativeCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceConservativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD conservative force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceConservativeCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceConservativeCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceConservativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD random force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD random force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceRandomCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceRandomCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceRandomCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD random force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD random force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceRandomCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceRandomCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceRandomCalculator: Unknown DPD type.");
        }
    }


    /**
     * Returns particle pair DPD dissipative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD dissipative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceDissipativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceDissipativeCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD dissipative force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD dissipative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceDissipativeCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceDissipativeCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }
    
    /**
     * Returns particle pair DPD potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairDpdPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdPotentialCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdPotentialCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD potential calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairDpdPotentialCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdPotentialCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdPotentialCalculator: Unknown DPD type.");
        }
    }


    /**
     * Returns particle-pair nearest-neighbor calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairNearestNeighborCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException 
    {
        return new ParticlePairNearestNeighborCalculator(
            this,
            aSimulationLogger, 
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength,
            aParallelizationInfo,
            aRandomNumberSeed
        );
    }
    // </editor-fold>
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- DPD parameters">
    /**
     * DPD cut-off length (in DPD units)
     * 
     * @return DPD cut-off length (in DPD units)
     */
    public float getDpdCutOffLength() {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                // DPD cut-off length of 1.0
                return 1.0f;
            default:
                throw new IllegalArgumentException("Factory.getDpdCutOffLength: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Harmonic bond calculators">
    /**
     * Returns harmonic bond conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @return Harmonic bond force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondForceCalculator getHarmonicBondForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        switch (this.bondType) {
            case HARMONIC:
                return new HarmonicBondForceConservativeCalculator(
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aParallelizationInfo
                );
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondForceConservativeCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond conservative force calculator
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @return Harmonic bond force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondForceCalculator getHarmonicBondForceConservativeCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        switch (this.bondType) {
            case HARMONIC:
                return new HarmonicBondForceConservativeCalculator(aHarmonicBondPropertyCalculator);
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondForceConservativeCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @return Harmonic bond potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondPropertyCalculator getHarmonicBondPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        switch (this.bondType) {
            case HARMONIC:
                return new HarmonicBondPotentialCalculator(
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aParallelizationInfo
                );
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondPotentialCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond potential calculator
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @return Harmonic bond potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondPropertyCalculator getHarmonicBondPotentialCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        switch (this.bondType) {
            case HARMONIC:
                return new HarmonicBondPotentialCalculator(aHarmonicBondPropertyCalculator);
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondPotentialCalculator: Unknown harmonic bond type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Electrostatics calculators">
    /**
     * Returns particle pair electrostatics conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair electrostatics force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairElectrostaticsForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocForceConservativeCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            case DPD:
                return new ParticlePairElectrostaticsDpdForceConservativeCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsForceConservativeCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics conservative force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair electrostatics force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairElectrostaticsForceConservativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocForceConservativeCalculator(aParticlePairInteractionCalculator);
            case DPD:
                return new ParticlePairElectrostaticsDpdForceConservativeCalculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsForceConservativeCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair electrostatics potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairElectrostaticsPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocPotentialCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            case DPD:
                return new ParticlePairElectrostaticsDpdPotentialCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsPotentialCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics potential calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair electrostatics potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairElectrostaticsPotentialCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocPotentialCalculator(aParticlePairInteractionCalculator);
            case DPD:
                return new ParticlePairElectrostaticsDpdPotentialCalculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsPotentialCalculator: Unknown electrostatics type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Time step calculators">
    /**
     * Time step calculator
     * 
     * @param aSimulationOutput Simulation output
     * @param aSimulationLogger Simulation logger
     * @param aParameters Parameters instance
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @param aParticlePositionPool Particle position pool (NOT allowed to be 
     * null)
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     * @return Time step calculator
     * @throws IllegalArgumentException Thrown if argument is invalid
     */
    public ITimeStepCalculator getTimeStepCalculator(
        IOutput aSimulationOutput,
        ILogger aSimulationLogger, 
        Parameters aParameters, 
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed,
        ParticlePositionPool aParticlePositionPool,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        boolean tmpIsCacheUsage = false;
        switch (this.integrationType) {
            case GWMVV:
                // Lambda parameter for Groot-Warren Modified Velocity-Verlet (GWMVV) integration
                float tmpLambda = (float) this.integrationParameters[0];
                return new GwmvvTimeStepCalculator(
                    this,
                    tmpLambda,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case SCMVV:
                // Number of self-constistent iterations for Self-consistent Modified Velocity-Verlet (SCMVV) integration
                int tmpSelfConsistentIterationNumber = (Integer) this.integrationParameters[0];
                // Flag for use of cache for Self-consistent Modified Velocity-Verlet (SCMVV) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[1];
                return new ScmvvTimeStepCalculator(
                    this,
                    tmpSelfConsistentIterationNumber,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case S1MVV:
                // Flag for use of cache for Shardlow S1 Modified Velocity-Verlet (S1MVV) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[0];
                return new S1mvvTimeStepCalculator(
                    this,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case PNHLN:
                // Mu parameter for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
                float tmpMu = (float) this.integrationParameters[0];
                // Flag for use of cache for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[1];
                return new PnhlnTimeStepCalculator(
                    this,
                    tmpMu,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            default:
                throw new IllegalArgumentException("Factory.getTimeStepCalculator: Unknown integration type.");
        }
    }
    // </editor-fold>
    // </editor-fold>
    
}

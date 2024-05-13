package ru.bmstu.metanapp;

import static ru.bmstu.metanapp.staticField.CNUMBR.*;
import static ru.bmstu.metanapp.staticField.CGED.*;
import static ru.bmstu.metanapp.staticField.CGEN.*;
import static ru.bmstu.metanapp.staticField.CVLT.*;
import static ru.bmstu.metanapp.staticField.CVLE.*;
import static ru.bmstu.metanapp.staticField.CVNT.*;
import static ru.bmstu.metanapp.staticField.XL.*;
import static ru.bmstu.metanapp.staticField.XK.*;
import static ru.bmstu.metanapp.staticField.CW.*;

public class HelloApplication {

    public static void main(String[] args) {
        //staticFields RK1 20
        trueInit(); // RK1 45
        double[][][] A = new double[46][22][16];
        int N1 = 46;
        int N2 = 22;
        int N3 = 16;
        RK2(N1, N2, N3, A);
    }

    public static void RK2(int N1, int N2, int N3, double[][][] A) {
        INM = IN - 1;
        JNM = JN - 1;
        NITER = 0;
        NVAR = 0;
        double S1;
        double G1;
        double DSS;
        double FSS;
        double FSKR;
        double FOS;
        double FA2;
        double APR;
        double AP2;
        double AP1;
        double NIT;
        double Z1;
        double Z2;
        double APQ;
        double ZAP1;
        double ZAP2;
        double ZAP3;
        double ZAP4;
        double FA1;
        double FZ;
        double APZ;
        double APO;
        double BETA = 0;
        double AKM1;
        double AKM;
        double FKR;
        double AMO2;
        double AMS1;
        double AKW;
        double AMO1;
        double AMG1;
        double ALPRIW;
        double VR;
        double ZGCP;
        double ZHP;
        double ZGCPE;
        double AMZ;
        double ZGCPT;
        double Z;
        do {
            NVAR++;
            ALFA += 0.2;
            double DELD = 0.5;
            if (NVAR != 1) {
                DELD += 1.E-3;
            }
            grid(N1, N2, N3, A);
            C = 3.14 * DOOK * DELOK;
            S1 = 0.785 * Math.pow(DCMB, 2);
            G1 = 3.14 * DOGOR * DELGOR;
            for (int I = 1; I <= IN; I++) {
                PKS[I - 1] = P;
                if (I < IC1) {
                    continue;
                }
                int JL = JC + IC - I;
                DSS = X2[JL - 1] * 2;
                FSS = 3.14 / 4. * Math.pow(DSS, 2);
                FSKR = 3.14 * Math.pow(X2[JC - 1], 2);
                FOS = FSS / FSKR;
                if (FOS <= 1) {
                    FOS = 1;
                }
                FA2 = 1. / FOS;
                APR = Math.pow((2. / (AKT + 1.)), (AKT / (AKT - 1.))) * P;
                AP2 = APR / P;
                AP1 = 1.;
                NIT = 0;
                while (true) {
                    NIT++;
                    Z1 = AP1;
                    Z2 = AP2;
                    APQ = (Z1 + Z2) / 2.;
                    ZAP1 = Math.pow((2. / (AKT + 1)), ((AKT + 1.) / (2. * (AKT - 1.))));
                    ZAP2 = Math.pow(((AKT - 1.) / 2.), 0.5);
                    ZAP3 = Math.pow((1. / APQ), (1. / AKT));
                    ZAP4 = Math.pow(1. - Math.pow(APQ, (AKT - 1.) / AKT), 0.5);
                    FA1 = 1. / (ZAP1 * ZAP2 * ZAP3 / ZAP4);
                    FZ = FA2 - FA1;
                    if (Math.abs(FZ) <= 1.E-6 || NIT > 1.E+4) {
                        break;
                    }
                    if (FA2 <= FA1) {
                        AP2 = APQ;
                    }
                    if (FA2 >= FA1) {
                        AP1 = APQ;
                    }
                }
                PKS[I - 1] = APQ * P;
                if (I > IC) {
                    PKS[I - 1] = PKS[IC - 1];
                }
            } // 80 RK2 84
            APZ = P - PKS[IC1 - 1];
            for (int I = 2; I <= IC1; I++) {
                PKS[I - 1] = P - APZ * X1[I - 1] / X1[IC1 - 1];
                APO = PKS[I - 1] / P;
            } // 80 RK2 89
            if (ALFA <= 1.) {
                BETA = 1794. + 491.252 * ALFA - 475.9317 * Math.pow(ALFA, 2);
            }
            if (ALFA > 1.) {
                BETA = 2154.733 - 402.252 * ALFA + 54.1698 * Math.pow(ALFA, 2);
            }
            if (ALFA < 0.6) {
                BETA = 1917. - 2500. * (0.6 - ALFA);
            }
            if (BETA < 600.) {
                BETA = 600.;
            }
            AKM1 = 70.0;
            DKR = 2.0 * X2[JC - 1];
            AKM = 4.0 * ALFA;
            FKR = Math.pow(DKR, 2) * Math.PI / 4.0;
            AMS = P * FKR / BETA;
            AMG = AMS / (1 + AKM);
            AMO = AMG * AKM;
            AMO2 = AMO;
            ROS = P / 77427.4;
            ROP = P / 154855.7;
            AMS1 = D;
            AKW = S1 / (C + S1);
            AMS1 = AKW * AMO;
            VINS1 = AMS1 / (ROS * Math.PI * Math.pow(R[JC - 1], 2));
            AMO1 = AMS1;
            AMG1 = 0.0;
            VINS = (AMO - AMO1) / ROS / (Math.PI * (Math.pow(R[JA2 - 1], 2) - Math.pow(R[JA1 - 1], 2)));
            VINP = (AMG - AMG1) / ROP / (Math.PI * (Math.pow(R[JA - 1], 2) - Math.pow(R[JZ1 - 1], 2)));
            DKS = X2[JB1 - 1] * 2.0;
            ALKS = X1[IC1 - 1] - X1[IC2 - 1];
            ALKD = X1[IC - 1];
            ALPRIW = (Math.pow(R[JA3 - 1], 2) * X1[IC2 - 1] + (Math.pow(R[JB1 - 1], 2) - Math.pow(R[JC - 1], 2)) * (X1[IC1 - 1] - X1[IC2 - 1]) + (Math.pow(R[JB1 - 1], 2) - Math.pow(R[JC - 1], 2)) * (X1[IC - 1] - X1[IC1 - 1]) / 3.0 + Math.pow(R[JC - 1], 2) * (X1[IC - 1] - X1[IC2 - 1])) / Math.pow(R[JC - 1], 2);
            TGY = (X2[JN - 1] - X2[JC - 1]) / (X1[IC - 1] - X1[IC1 - 1]);
            DKAR = X2[JA3 - 1] * 2.0;
            DELKAR = X1[IC2 - 1];
            DOOK = X2[JA1 - 1] * 2.0;
            DELOK = X2[JA2 - 1] - X2[JA1 - 1];
            DOGOR = X2[JZ1 - 1] * 2.0;
            DELGOR = X2[JA - 1] - X2[JZ1 - 1];
            DCMB = X2[JZ - 1] * 2.0;
            init(N1, N2, N3, A);
            VR = A[1 - 1][JZ1 + 1 - 1][NMU2 - 1] / R[JZ1 + 1 - 1];
            write7("RASHET N " + NVAR);
            write7(String.format("P=%.4e ALFA=%.2f KM=%.2f BETA=%.1f Ms=%.4e Mok=%.4e Mgor=%.4e Vgor=%.1f Vok=%.1f ROg=%.1f ROok=%.1f Vokr=%.1f Vcmw=%.1f Mcmw=%.4e%n",
                    P, ALFA, AKM, BETA, AMS, AMO, AMG, VINP, VINS, ROP, ROS, VR, VINS1, AMS1));

            write7(String.format("JZ=%d JZ1=%d JA=%d JA1=%d JA2=%d JA3=%d JC=%d JN=%d IC=%d IC1=%d IC2=%d IN=%d%n",
                    JZ, JZ1, JA, JA1, JA2, JA3, JC, JN, IC, IC1, IC2, IN));

            write7(String.format("Dks=%.4e Lks=%.4e Lkd=%.4e Dkr=%.4e tgY=%.1f Dprk=%.4e Lprk=%.4e Dook=%.4e DELOK=%.4e DOGOR=%.4e DELGOR=%.4e Dcmw=%.4e Lpriw=%.4e Fod=%.4e Fcmb=%.4e Fgor=%.4e Kwosp=%.4e%n",
                    DKS, ALKS, ALKD, DKR, TGY, DKAR, DELKAR, DOOK, DELOK, DOGOR, DELGOR, DCMB, ALPRIW, C, S1, G1, AKW));
            NITER = 0;
            do {
                NITER++;
                viscos(N1, N2, N3, A);
                veldis(N1, N2, N3, A);
                convec(N1, N2, N3, A);
                eqn(N1, N2, N3, A);
            } while (NITER >= NMAX);
            write7("NITER=" + NITER);
            int JC17 = JC - 1;
            if (ALFA <= 1) {
                //goto 28
                ZGCP = 0;
                for (int J = 2; J <= JC17; J++) {
                    ZGCP = ZGCP + A[IN - 1][J - 1][NMFU1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2;
                }
                ZGCP = ZGCP + A[IN - 1][1 - 1][NMFU1 - 1] * (A[IN - 1][2 - 1][NF - 1] - A[IN - 1][1 - 1][NF - 1]) / 2. + A[IN - 1][JC - 1][NMFU1 - 1] * (A[IN - 1][JC - 1][NF - 1] - A[IN - 1][JC17 - 1][NF - 1]) / 2.;
                AMZ = AMG - AMO2 / 4; //RK2 231
                ZGCP = (ZGCP - AMZ / 6.28) / ((AMG - AMZ) / 6.28);
                if (Math.abs(ZGCP) > 1) {
                    ZGCP = 1;
                }
                ZGCPT = Math.sqrt((1 - Math.abs(ZGCP)));
                write7(String.format("************* KOEF Fkomp (OK) = %.5f ***********%n", ZGCPT));
            } else {
                ZGCP = 0;
                ZHP = 0;
                for (int J = 2; J <= JC17; J++) {
                    ZHP = ZHP + A[IN - 1][J - 1][NMU1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2;
                    ZGCP = ZGCP + A[IN - 1][J - 1][NMFU1 - 1] * (A[IN - 1][J + 1 - 1][NF - 1] - A[IN - 1][J - 1 - 1][NF - 1]) / 2;
                }
                ZGCP += A[IN - 1][1 - 1][NMFU1 - 1] * (A[IN - 1][2 - 1][NF - 1] - A[IN - 1][1 - 1][NF - 1]) / 2.0 +
                        A[IN - 1][JC - 1][NMFU1 - 1] * (A[IN - 1][JC - 1][NF - 1] - A[IN - 1][JC17 - 1][NF - 1]) / 2.0;
                ZGCP /= ((A[1 - 1][JA - 1][NF - 1] - A[1 - 1][JZ1 - 1][NF - 1]) + 0.013833 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]));
                ZHP /= (A[1 - 1][JC - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);
                if (Math.abs(ZGCP) > 1) {
                    ZGCP = 1;
                }
                ZGCPE = Math.sqrt(1 - Math.abs(ZGCP)); //RK2 218
                write7(String.format("************* KOEF Fkomp (GOR) = %.5f **********%n", ZGCPE));
            }
            Z = 0;
            double ZO;
            double ZG;
            double ALFAJ;
            double BETAJ;
            double ZF;
            for (int J = 2; J <= JC; J++) {
                int IL = INM;
                ZO = A[IL - 1][J - 1 - 1][NMOX1 - 1] + 4.27 * A[IL - 1][J - 1 - 1][NMPR1 - 1] + 3.35 * A[IL - 1][J - 1 - 1][NMOX2 - 1] +
                        376. * A[IL - 1][J - 1 - 1][NMDF - 1] + 2.785 * A[IL - 1][J - 1 - 1][NMPR2 - 1] + 2.765 * A[IL - 1][J - 1 - 1][NMFU2 - 1];
                ZG = A[IL - 1][J - 1 - 1][NMFU1 - 1] + 1.09 * A[IL - 1][J - 1 - 1][NMPR1 - 1] + 0.857 * A[IL - 1][J - 1 - 1][NMOX2 - 1] +
                        96. * A[IL - 1][J - 1 - 1][NMDF - 1] + 0.711 * A[IL - 1][J - 1 - 1][NMPR2 - 1] + 0.706 * A[IL - 1][J - 1 - 1][NMFU2 - 1];
                ALFAJ = ZO / ZG / 3.9166;
                if (ALFAJ <= 1.0) {
                    BETAJ = 1794. + 491.252 * ALFAJ - 475.9317 * Math.pow(ALFAJ, 2);
                } else {
                    BETAJ = 2154. - 402.252 * ALFAJ + 54.1698 * Math.pow(ALFAJ, 2);
                }
                if (ALFAJ < 0.6)
                    BETAJ = 1917. - 2500. * (0.6 - ALFAJ);
                if (BETAJ < 600.0)
                    BETAJ = 600.0;
                ZF = (A[IL - 1][J - 1][NF - 1] - A[IL - 1][J - 1 - 1][NF - 1]) / (A[IL - 1][JC - 1][NF - 1] - A[IL - 1][1 - 1][NF - 1]);
                Z = Z + BETAJ * ZF;
            }

            double FICM = Z / BETA;
            System.out.printf("%25s FIbeta=%11.4E", "", FICM);
            print1(N1, N2, N3, A, IV);
            double FIOP = A[1 - 1][JA - 1][NF - 1] - A[1 - 1][JZ1 - 1][NF - 1] + 0.013833 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);
            double FIOS = A[IC2 - 1][JA2 - 1][NF - 1] - A[IC2 - 1][JA1 - 1][NF - 1] + 0.986167 * (A[1 - 1][JZ - 1][NF - 1] - A[1 - 1][1 - 1][NF - 1]);
            FIOP = FIOP * 3.1416 * 2;
            FIOS = FIOS * 3.1416 * 2;
            double DOP = FIOS / FIOP;
        } while (NVAR != NVARM);
    }

    private static void print1(int N1, int N2, int N3, double[][][] A, int IV) {
        String format101 = "%2s%11.3f%11.3f%11.3f%11.3f%11.3f%11.3f%n";
        String format102 = "%1s     %2s         %3s         %4s         %5s         %6s    X1(I)      %7s      X2(J)%n";

        write9(String.format(format102, "", "2", "3", "4", "5", "6", "9", "10", "11", "12", "13"));
        for (int J = 0; J < JN; J++) {
            for (int I = 0; I < IN; I++) {
                write9(String.format(format101, "", A[I][J][2], A[I][J][3], A[I][J][4], A[I][J][5], A[I][J][6], X1[I], X2[J]));
            }
        }
        // Записываем вторую строку данных
        write9(String.format(format102, "", "2", "3", "4", "5", "6", "9", "10", "11", "12", "13"));
        for (int J = 0; J < JN; J++) {
            for (int I = 0; I < IN; I++) {
                write9(String.format(format101, "", A[I][J][9], A[I][J][10], A[I][J][11], A[I][J][12], A[I][J][13], X1[I], X2[J]));
            }
        }

    }

    private static void eqn(int N1, int N2, int N3, double[][][] A) {
        double SOURCE = 0;
        double RSO = 0;
        double Z = 0;
        double BBE;
        double ADNM;
        double BBW;
        double BBS;
        double BBN;
        double AAE;
        double AAW;
        double AAN;
        double AAS;
        double ANUM;
        double RISO;
        double ROP3;
        double RS = 0;
        if (NITER <= NPRINT && NVAR == 1) {
            //goto 177
        }
        RL = 0.1;
        for (int J = 2; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) {
                if (
                        (J == JA && I == 2) ||
                                (J == JA1 && I == (IC2 + 1)) ||
                                (J == JA2 && I == (IC2 + 1)) ||
                                (J == JA3 && I == (IC2 + 1)) ||
                                (J == JZ && I == 2) ||
                                (J == JZ1 && I == 2)
                ) {
                    Z = A[I - 1][J - 1][NW - 1]; // RK3 55
                    RISO = 1.0 / (R[J - 1] * R[J - 1]);
                    ROP3 = A[I - 1][J - 1][NRO3 - 1];
                    BBE = 4.0 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                    BBW = 4.0 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                    BBN = 16.0 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / ((R[J + 1 - 1] + R[J - 1]) * (R[J + 1 - 1] + R[J - 1])) * BN[J - 1];
                    BBS = 16.0 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / ((R[J - 1 - 1] + R[J - 1]) * (R[J - 1 - 1] + R[J - 1])) * BS[J - 1];
                    ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                    ADNM = BBE + BBW + BBN + BBS;
                    A[I - 1][J - 1][NW - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                    A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    //goto 111
                } else if (JB1 == JN) {
                    Z = A[I - 1][J - 1][NW - 1];
                    SOURCE = sorce(N1, N2, N3, A, SOURCE, I, J, NW);
                    RSO = R[J - 1] * R[J - 1];
                    BBE = 2.0 * RSO * BE[I - 1];
                    BBW = 2.0 * RSO * BW[I - 1];
                    BBS = (R[J - 1 - 1] * R[J - 1 - 1] + RSO) * BS[J - 1];
                    BBN = (R[J + 1 - 1] * R[J + 1 - 1] + RSO) * BN[J - 1];
                    AAE = RSO * AE[I - 1][J - 1];
                    AAW = RSO * AW[I - 1][J - 1];
                    AAN = RSO * AN[I - 1][J - 1];
                    AAS = RSO * AS[I - 1][J - 1];
                    ANUM = (AAE + A[I + 1 - 1][J - 1][NMU - 1] * BBE) * A[I + 1 - 1][J - 1][NW - 1] +
                            (AAW + A[I - 1 - 1][J - 1][NMU - 1] * BBW) * A[I - 1 - 1][J - 1][NW - 1] +
                            (AAN + A[I - 1][J + 1 - 1][NMU - 1] * BBN) * A[I - 1][J + 1 - 1][NW - 1] +
                            (AAS + A[I - 1][J - 1 - 1][NMU - 1] * BBS) * A[I - 1][J - 1 - 1][NW - 1] + SOURCE;
                    ADNM = AAE + AAW + AAN + AAS + A[I - 1][J - 1][NMU - 1] * (BBE + BBW + BBN + BBS);
                    if (ADNM != 0) {
                        //goto 112
                        A[I - 1][J - 1][NW - 1] = ANUM / ADNM;
                        A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    }
                    //goto 77
                } else if ((J == (JB1 - 1)) && (I == IC2)) {
                    Z = A[I - 1][J - 1][NW - 1];
                    RISO = 1 / R[J - 1] / R[J - 1];
                    ROP3 = A[I - 1][J - 1][NRO3 - 1];
                    BBE = 4 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                    BBW = 4 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                    BBN = 16 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J + 1 - 1] + R[J - 1], 2)) * BN[J - 1];
                    BBS = 16 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J - 1 - 1] + R[J - 1], 2)) * BS[J - 1];
                    ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                    ADNM = BBE + BBW * BBN * BBS;
                    A[I - 1][J - 1][NF - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                    A[I - 1][J - 1][NF - 1] = Z + (A[I - 1][J - 1][NF - 1] - Z) * RL;
                    //goto 113
                } else {
                    Z = A[I - 1][J - 1][NW - 1];
                    SOURCE = sorce(N1, N2, N3, A, SOURCE, I, J, NW);
                    RSO = R[J - 1] * R[J - 1];
                    BBE = 2.0 * RSO * BE[I - 1];
                    BBW = 2.0 * RSO * BW[I - 1];
                    BBS = (R[J - 1 - 1] * R[J - 1 - 1] + RSO) * BS[J - 1];
                    BBN = (R[J + 1 - 1] * R[J + 1 - 1] + RSO) * BN[J - 1];
                    AAE = RSO * AE[I - 1][J - 1];
                    AAW = RSO * AW[I - 1][J - 1];
                    AAN = RSO * AN[I - 1][J - 1];
                    AAS = RSO * AS[I - 1][J - 1];
                    ANUM = (AAE + A[I + 1 - 1][J - 1][NMU - 1] * BBE) * A[I + 1 - 1][J - 1][NW - 1] +
                            (AAW + A[I - 1 - 1][J - 1][NMU - 1] * BBW) * A[I - 1 - 1][J - 1][NW - 1] +
                            (AAN + A[I - 1][J + 1 - 1][NMU - 1] * BBN) * A[I - 1][J + 1 - 1][NW - 1] +
                            (AAS + A[I - 1][J - 1 - 1][NMU - 1] * BBS) * A[I - 1][J - 1 - 1][NW - 1] + SOURCE;
                    ADNM = AAE + AAW + AAN + AAS + A[I - 1][J - 1][NMU - 1] * (BBE + BBW + BBN + BBS);
                    if (ADNM != 0) {
                        //goto 112
                        A[I - 1][J - 1][NW - 1] = ANUM / ADNM;
                        A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                    }
                }
                /*// 77
                Z = A[I - 1][J - 1][NW - 1];
                sorce(N1, N2, N3, A, SOURCE, I, J, NW);
                RSO = R[J - 1] * R[J - 1];
                BBE = 2.0 * RSO * BE[I - 1];
                BBW = 2.0 * RSO * BW[I - 1];
                BBS = (R[J - 1 - 1] * R[J - 1 - 1] + RSO) * BS[J - 1];
                BBN = (R[J + 1 - 1] * R[J + 1 - 1] + RSO) * BN[J - 1];
                AAE = RSO * AE[I - 1][J - 1];
                AAW = RSO * AW[I - 1][J - 1];
                AAN = RSO * AN[I - 1][J - 1];
                AAS = RSO * AS[I - 1][J - 1];
                ANUM = (AAE + A[I + 1 - 1][J - 1][NMU - 1] * BBE) * A[I + 1 - 1][J - 1][NW - 1] +
                        (AAW + A[I - 1 - 1][J - 1][NMU - 1] * BBW) * A[I - 1 - 1][J - 1][NW - 1] +
                        (AAN + A[I - 1][J + 1 - 1][NMU - 1] * BBN) * A[I - 1][J + 1 - 1][NW - 1] +
                        (AAS + A[I - 1][J - 1 - 1][NMU - 1] * BBS) * A[I - 1][J - 1 - 1][NW - 1] + SOURCE;
                ADNM = AAE + AAW + AAN + AAS + A[I - 1][J - 1][NMU - 1] * (BBE + BBW + BBN + BBS);
                if (ADNM == 0) {
                    //goto 112
                }
                A[I - 1][J - 1][NW - 1] = ANUM / ADNM;
                A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                //112
                //goto 114
                //111
                Z = A[I - 1][J - 1][NW - 1]; // RK3 55
                RISO = 1.0 / (R[J - 1] * R[J - 1]);
                ROP3 = A[I - 1][J - 1][NRO3 - 1];
                BBE = 4.0 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                BBW = 4.0 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                BBN = 16.0 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / ((R[J + 1 - 1] + R[J - 1]) * (R[J + 1 - 1] + R[J - 1])) * BN[J - 1];
                BBS = 16.0 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / ((R[J - 1 - 1] + R[J - 1]) * (R[J - 1 - 1] + R[J - 1])) * BS[J - 1];
                ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                ADNM = BBE + BBW + BBN + BBS;
                A[I - 1][J - 1][NW - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z) * RL;
                //goto 114
                //113 RK3 68
                Z = A[I - 1][J - 1][NW - 1];
                RISO = 1 / R[J - 1] / R[J - 1];
                ROP3 = A[I - 1][J - 1][NRO3 - 1];
                BBE = 4 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                BBW = 4 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                BBN = 16 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J + 1 - 1] + R[J - 1], 2)) * BN[J - 1];
                BBS = 16 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J - 1 - 1] + R[J - 1], 2)) * BS[J - 1];
                ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1];
                ADNM = BBE + BBW * BBN * BBS;
                A[I - 1][J - 1][NF - 1] = A[I - 1][J - 1][NF - 1] * ADNM - ANUM;
                A[I - 1][J - 1][NF - 1] = Z + (A[I - 1][J - 1][NF - 1] - Z) * RL;*/
                //114

                if (Math.abs(A[I - 1][J - 1][NW - 1]) >= 1.E-20) {
                    RS = 1 - Z / A[I - 1][J - 1][NW - 1];
                    //goto 15
                } else {
                    if (Math.abs(Z) >= 1.E-20) {
                        RS = 1 - A[I - 1][J - 1][NW - 1] / Z;
                    }
                }
                /*if (Math.abs(Z) >= 1.E-20) {
                    RS = 1 - A[I - 1][J - 1][NW - 1] / Z;
                }
                //goto 16
                RS = 1 - Z / A[I - 1][J - 1][NW - 1]; //15
                //16*/
                if (Math.abs(A[I - 1][J - 1][NW - 1]) < 1.E-20 && Math.abs(Z) < 1.E-20) {
                    RS = 0;
                }
                if (Math.abs(RS) > Math.abs(RSDU[NW - 1])) {
                    RSDU[NW - 1] = RS;
                }
                if (NITER >= 20) {
                    A[I - 1][J - 1][NW - 1] = Z + (A[I - 1][J - 1][NW - 1] - Z);
                }
                if (A[I - 1][J - 1][NW - 1] >= 1.E+14) {
                    A[I - 1][J - 1][NW - 1] = 1.E+14;
                }
                if (A[I - 1][J - 1][NW - 1] <= -1.E+14) {
                    A[I - 1][J - 1][NW - 1] = -1.E+14;
                }
                //11
                //177 RK3 93
            }
        }
        for (int J = 2; J <= JNM; J++) { // to 21
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) {
                Z = A[I - 1][J - 1][NF - 1];
                if ((J == JA && I == 2) || (J == JA1 && I == IC2 + 1) || (J == JA2 && I == IC2 + 1) || (J == JA3 && I == IC2 + 1) || (J == JZ && I == 2) || (J == JZ1 && I == 2)) {
                    A[I - 1][J - 1][NF - 1] = A[I - 1 - 1][J - 1][NF - 1];
                    //goto 121
                } else if (JB1 != JN && J == (JB1 - 1) && I == IC2) {
                    A[I - 1][J - 1][NF - 1] = A[I - 1][J + 1 - 1][NF - 1];
                    //goto 123
                } else {
                    SOURCE = sorce(N1, N2, N3, A, SOURCE, I, J, NF);
                    if (NITER <= NPRINT) {
                        SOURCE = 0;
                    }
                    RISO = 1 / R[J - 1] / R[J - 1];
                    ROP3 = A[I - 1][J - 1][NRO3 - 1];
                    BBE = 4 / (A[I + 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BE[I - 1];
                    BBW = 4 / (A[I - 1 - 1][J - 1][NRO3 - 1] + ROP3) * RISO * BW[I - 1];
                    BBN = 16 / (A[I - 1][J + 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J + 1 - 1] + R[J - 1], 2)) * BN[J - 1];
                    BBS = 16 / (A[I - 1][J - 1 - 1][NRO3 - 1] + ROP3) / (Math.pow(R[J - 1 - 1] + R[J - 1], 2)) * BS[J - 1];
                    ANUM = BBE * A[I + 1 - 1][J - 1][NF - 1] + BBW * A[I - 1 - 1][J - 1][NF - 1] + BBN * A[I - 1][J + 1 - 1][NF - 1] + BBS * A[I - 1][J - 1 - 1][NF - 1] + SOURCE;
                    ADNM = BBE + BBW + BBN + BBS;
                    if (ADNM == 0) {
                        continue;
                    }
                    A[I - 1][J - 1][NF - 1] = ANUM / ADNM;
                }

                if (Math.abs(A[I - 1][J - 1][NF - 1]) >= 1.E-20) {
                    RS = 1 - Z / A[I - 1][J - 1][NF - 1];
                    //goto 25
                } else {
                    if (Math.abs(Z) >= 1.E-20) {
                        RS = 1 - A[I - 1][J - 1][NF - 1] / Z;
                    }
                }
                if (Math.abs(A[I - 1][J - 1][NF - 1]) < 1.E-20 && Math.abs(Z) < 1.E-20) {
                    RS = 0;
                }
                if (Math.abs(RS) > Math.abs(RSDU[NF - 1])) {
                    RSDU[NF - 1] = RS;
                }
                if (A[I - 1][J - 1][NF - 1] >= 1.E+25) {
                    A[I - 1][J - 1][NF - 1] = 1.E+25;
                }
                if (A[I - 1][J - 1][NF - 1] <= -1.E+25) {
                    A[I - 1][J - 1][NF - 1] = -1.E+25;
                }
            }
        } // 21 RK3 137
        for (int J = 3; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) {
                BBE = (A[I + 1 - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2) + A[I - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2)) * BE[I - 1];
                BBW = (A[I - 1 - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2) + A[I - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2)) * BW[I - 1];
                BBN = (A[I - 1][J + 1 - 1][NMU - 1] * Math.pow(R[J - 1], 2) + A[I - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2)) * BN[I - 1];
                BBS = (A[I - 1][J - 1 - 1][NMU - 1] * Math.pow(R[J - 1], 2) + A[I - 1][J - 1][NMU - 1] * Math.pow(R[J - 1], 2)) * BS[I - 1];
                ANUM = (AE[I - 1][J - 1] + BBE / Math.pow(R[J - 1], 2)) * A[I + 1 - 1][J - 1][NMU2 - 1] +
                        (AW[I - 1][J - 1] + BBW / Math.pow(R[J - 1], 2)) * A[I - 1 - 1][J - 1][NMU2 - 1] +
                        (AN[I - 1][J - 1] + BBN / Math.pow(R[J + 1 - 1], 2)) * A[I - 1][J + 1 - 1][NMU2 - 1] +
                        (AS[I - 1][J - 1] + BBS / Math.pow(R[J - 1 - 1], 2)) * A[I - 1][J - 1 - 1][NMU2 - 1];
                ADNM = AE[I - 1][J - 1] + AW[I - 1][J - 1] + AN[I - 1][J - 1] + AS[I - 1][J - 1] +
                        BBE / Math.pow(R[J - 1], 2) + BBW / Math.pow(R[J - 1], 2) + BBN / Math.pow(R[J + 1 - 1], 2) +
                        BBS / Math.pow(R[J - 1 - 1], 2);
                if (ADNM == 0) {
                    continue;
                }
                Z = A[I - 1][J - 1][NMU2 - 1];
                A[I - 1][J - 1][NMU2 - 1] = ANUM / ADNM;
            }
        }//RK3 156

        for (int I = 2; I <= INM; I++) {
            A[I - 1][2 - 1][NMU2 - 1] = A[I - 1][3 - 1][NMU2 - 1] * R[2 - 1] / R[3 - 1];
        }

        for (int II = 3; II <= 10; II++) {
            int K = 13 - II;
            for (int J = 2; J <= JNM; J++) {
                int IL = IMIN[J - 1];
                int IM = IMAX[J - 1];
                if (J == JA3) {
                    IL = IC2 + 1;
                }
                if (J == JC) {
                    IM = IC - 1;
                }
                for (int I = IL; I <= IM; I++) {
                    ADNM = AE1[I - 1][J - 1] + AW1[I - 1][J - 1] + AN1[I - 1][J - 1] + AS1[I - 1][J - 1];
                    Z = A[I - 1][J - 1][K - 1]; // RK3 169
                    ANUM = AE1[I - 1][J - 1] * A[I + 1 - 1][J - 1][K - 1] + AW1[I - 1][J - 1] * A[I - 1 - 1][J - 1][K - 1] + AN1[I - 1][J - 1] * A[I - 1][J + 1 - 1][K - 1] + AS1[I - 1][J - 1] * A[I - 1][J - 1 - 1][K - 1];
                    SOURCE = sorce(N1, N2, N3, A, SOURCE, I, J, K);
                    ANUM = ANUM + SOURCE;
                    if (ADNM == 0) {
                        continue;
                        //goto 51
                    }
                    A[I - 1][J - 1][K - 1] = ANUM / ADNM;
                    if (K != NMU1) {
                        continue;
                        //goto 51
                    }
                    if (Math.abs(A[I - 1][J - 1][K - 1]) >= 1.E-20) {
                        RS = 1 - Z / A[I - 1][J - 1][K - 1]; // 65
                        //goto 65
                    } else {
                        if (Math.abs(Z) >= 1.E-20) {
                            RS = 1 - A[I - 1][J - 1][K - 1] / Z;
                        }
                    }
                    if (Math.abs(A[I - 1][J - 1][K - 1]) < 1.E-20 && Math.abs(Z) < 1.E-20) {
                        RS = 0;
                    }
                    if (Math.abs(RS) > Math.abs(RSDU[K - 1])) {
                        RSDU[K - 1] = RS;
                    }
                } //51
            }//51
        }//61
        double CM;
        double ZNK;
        for (int J = 2; J <= JNM; J++) { //to 56
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) { //to 56
                FKS[I - 1][J - 1] = A[I - 1][J - 1][NMFU1 - 1];
                TKS[I - 1][J - 1] = A[I - 1][J - 1][NMOX1 - 1];
                CM = A[I - 1][J - 1][NMFU1 - 1] + A[I - 1][J - 1][NMFU2 - 1] + A[I - 1][J - 1][NMOX1 - 1] + A[I - 1][J - 1][NMOX2 - 1] + A[I - 1][J - 1][NMPR1 - 1] + A[I - 1][J - 1][NMPR2 - 1] + A[I - 1][J - 1][NMDF - 1];
                if (Math.abs(CM) <= 1.E-20) {
                    continue;
                }
                A[I - 1][J - 1][NMFU1 - 1] = A[I - 1][J - 1][NMFU1 - 1] / CM;
                A[I - 1][J - 1][NMFU2 - 1] = A[I - 1][J - 1][NMFU2 - 1] / CM;
                A[I - 1][J - 1][NMOX1 - 1] = A[I - 1][J - 1][NMOX1 - 1] / CM;
                A[I - 1][J - 1][NMOX2 - 1] = A[I - 1][J - 1][NMOX2 - 1] / CM;
                A[I - 1][J - 1][NMPR1 - 1] = A[I - 1][J - 1][NMPR1 - 1] / CM;
                A[I - 1][J - 1][NMPR2 - 1] = A[I - 1][J - 1][NMPR2 - 1] / CM;
                A[I - 1][J - 1][NMDF - 1] = A[I - 1][J - 1][NMDF - 1] / CM;
            }
        }
        for (int J = 2; J <= JNM; J++) { // to 82
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1]; // RK3 207
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) { // to 82
                boolean goto331 = false;
                boolean goto33 = false;
                boolean goto304 = false;
                if (NITER != 18) {
                    goto331 = true;
                }
                if (I >= 3 && I <= 11 && J > JZ && J < JA1 && !goto331) {
                    goto33 = true;
                }
                if (A[I - 1][J - 1][NT - 1] < ROM1 && goto331 && !goto33) { //331
                    goto304 = true;
                }
                //33
                if (goto33 && !goto304) {
                    if (A[I - 1][J - 1][NMFU1 - 1] <= 1.E-20) {
                        continue;
                    }
                    ZNK = A[I - 1][J - 1][NMOX1 - 1] / A[I - 1][J - 1][NMFU1 - 1];
                    if (ZNK > 3.9166) {
                        A[I - 1][J - 1][NMPR1 - 1] = A[I - 1][J - 1][NMPR1 - 1] + (8. * 44. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        A[I - 1][J - 1][NMOX2 - 1] = A[I - 1][J - 1][NMOX2 - 1] + (28. * 16. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        A[I - 1][J - 1][NMDF - 1] = A[I - 1][J - 1][NMDF - 1] + (4. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        A[I - 1][J - 1][NMPR2 - 1] = A[I - 1][J - 1][NMPR2 - 1] + (30. * 18. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        A[I - 1][J - 1][NMFU2 - 1] = A[I - 1][J - 1][NMFU2 - 1] + (32. * 17. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        A[I - 1][J - 1][NMOX1 - 1] = A[I - 1][J - 1][NMOX1 - 1] - (47. * 32. / 24. / 16.) * A[I - 1][J - 1][NMFU1 - 1];
                        if (A[I - 1][J - 1][NMOX1 - 1] < 1.E-20) {
                            A[I - 1][J - 1][NMOX1 - 1] = 1.E-20;
                        }
                        A[I - 1][J - 1][NMFU1 - 1] = 1.E-20;
                    } else {
                        A[I - 1][J - 1][NMPR1 - 1] = A[I - 1][J - 1][NMPR1 - 1] + (8. * 44. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        A[I - 1][J - 1][NMOX2 - 1] = A[I - 1][J - 1][NMOX2 - 1] + (16. * 28. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        A[I - 1][J - 1][NMDF - 1] = A[I - 1][J - 1][NMDF - 1] + (4. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        A[I - 1][J - 1][NMPR2 - 1] = A[I - 1][J - 1][NMPR2 - 1] + (30. * 18. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        A[I - 1][J - 1][NMFU2 - 1] = A[I - 1][J - 1][NMFU2 - 1] + (32. * 17. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        A[I - 1][J - 1][NMFU1 - 1] = A[I - 1][J - 1][NMFU1 - 1] - (24. * 16. / 32. / 47.) * A[I - 1][J - 1][NMOX1 - 1];
                        if (A[I - 1][J - 1][NMFU1 - 1] < 1.E-20) {
                            A[I - 1][J - 1][NMFU1 - 1] = 1.E-20;
                        }
                        A[I - 1][J - 1][NMOX1 - 1] = 1.E-20;
                    }
                }
                //304 RK3 237
                CM = A[I - 1][J - 1][NMFU1 - 1] + A[I - 1][J - 1][NMFU2 - 1] + A[I - 1][J - 1][NMOX1 - 1] + A[I - 1][J - 1][NMOX2 - 1] + A[I - 1][J - 1][NMPR1 - 1] + A[I - 1][J - 1][NMPR2 - 1] + A[I - 1][J - 1][NMDF - 1];
                if (Math.abs(CM) <= 1.E-20) {
                    continue;
                }
                A[I - 1][J - 1][NMFU1 - 1] = A[I - 1][J - 1][NMFU1 - 1] / CM;
                A[I - 1][J - 1][NMFU2 - 1] = A[I - 1][J - 1][NMFU2 - 1] / CM;
                A[I - 1][J - 1][NMOX1 - 1] = A[I - 1][J - 1][NMOX1 - 1] / CM;
                A[I - 1][J - 1][NMOX2 - 1] = A[I - 1][J - 1][NMOX2 - 1] / CM;
                A[I - 1][J - 1][NMPR1 - 1] = A[I - 1][J - 1][NMPR1 - 1] / CM;
                A[I - 1][J - 1][NMPR2 - 1] = A[I - 1][J - 1][NMPR2 - 1] / CM;
                A[I - 1][J - 1][NMDF - 1] = A[I - 1][J - 1][NMDF - 1] / CM;
                if (Math.abs(FKS[I - 1][J - 1]) < 1.E-22) {
                    RS = 1;
                } else {
                    RS = 1 - A[I - 1][J - 1][NMFU1 - 1] / FKS[I - 1][J - 1];
                }
                if (Math.abs(FKS[I - 1][J - 1]) < 1.E-22 && Math.abs(A[I - 1][J - 1][NMFU1 - 1]) < 1.E-22) {
                    RS = 0;
                }
                if (Math.abs(RS) >= Math.abs(RSDU[NMFU1 - 1])) {
                    RSDU[NMFU1 - 1] = Math.abs(RS);
                }
                if (Math.abs(TKS[I - 1][J - 1]) < 1.E-22) {
                    RS = 1; //393
                } else {
                    RS = 1 - A[I - 1][J - 1][NMOX1 - 1] / TKS[I - 1][J - 1];
                }
                if (Math.abs(TKS[I - 1][J - 1]) < 1.E-22 && Math.abs(A[I - 1][J - 1][NMOX1 - 1]) < 1.E-22) {
                    RS = 0;
                }
                if (Math.abs(RS) >= Math.abs(RSDU[NMOX1 - 1])) {
                    RSDU[NMOX1 - 1] = Math.abs(RS);
                }
                FKS[I - 1][J - 1] = 0;
                TKS[I - 1][J - 1] = 0;
            }
        }
        bound(N1, N2, N3, A);
    }

    private static void bound(int N1, int N2, int N3, double[][][] A) {
        double WX;
        double WY;
        double WZ;
        double Z;
        int JA17 = JA + 1;
        int JA18 = JA1 - 1;
        RL = 1.;
        for (int J = JA17; J <= JN; J++) {
            int IL = IMIN[J-1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            for (int K = 3; K <= 10; K++) {
                if (J > JA18 && J <= JA2) {
                    continue; // Пропускаем итерацию, если условие выполняется
                }
                A[IL - 1-1][J-1][K-1] = A[IL-1][J-1][K-1];
            }
        }
        for (int J = JA; J <= JNM; J++) {
            int IL = IMIN[J-1];
            if (J == JA3) {
                IL = IC2 + 1;
            }
            if (J > JA1 && J < JA2) {
                continue;
            }

            WX = (A[IL-1][J-1][NF-1] - A[IL - 1-1][J-1][NF-1]) / Math.pow((R[J-1] * X1[IL-1]), 2);
            WY = A[IL-1][J-1][NW-1] * ((A[IL - 1-1][J-1][NRO3-1]) / 6.0 + (A[IL-1][J-1][NRO3-1] - A[IL - 1-1][J-1][NRO3-1]) / 8.0);
            WZ = A[IL - 1-1][J-1][NRO3-1] / 3.0 + 5.0 * (A[IL-1][J-1][NRO3-1] - A[IL - 1-1][J-1][NRO3-1]) / 24.0;
            A[IL - 1-1][J-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (A[IL - 1-1][J-1][NW-1] >= 1.0E+14) {
                A[IL - 1-1][J-1][NW-1] = 1.0E+14;
            }
            if (A[IL - 1-1][J-1][NW-1] <= -1.0E+14) {
                A[IL - 1-1][J-1][NW-1] = -1.0E+14;
            }
            if (NITER <= NPRINT) {
                A[IL - 1-1][J-1][NW-1] = 0.0;
            }
        }
        for (int J = JZ; J <= JZ1; J++) {
            WX = (A[2-1][J-1][NF-1] - A[1-1][J-1][NF-1]) / Math.pow((R[J-1] * X1[2-1]), 2);
            WY = A[2-1][J-1][NW-1] * ((A[1-1][J-1][NRO3-1]) / 6.0 + (A[2-1][J-1][NRO3-1] - A[1-1][J-1][NRO3-1]) / 8.0);
            WZ = A[1-1][J-1][NRO3-1] / 3.0 + 5.0 * (A[2-1][J-1][NRO3-1] - A[1-1][J-1][NRO3-1]) / 24.0;
            A[1-1][J-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (A[1-1][J-1][NW-1] >= 1.0E+14) {
                A[1-1][J-1][NW-1] = 1.0E+14;
            }
            if (A[1-1][J-1][NW-1] <= -1.0E+14) {
                A[1-1][J-1][NW-1] = -1.0E+14;
            }
            if (NITER <= NPRINT && NVAR == 1) {
                A[1-1][J-1][NW-1] = 0.0;
            }
        }
        int JZ17 = JZ + 1;
        int JZ18 = JZ1 - 1;
        for (int J = JZ17; J <= JZ18; J++) {
            for (int K = 3; K <= 10; K++) {
                A[1-1][J-1][K-1] = A[2-1][J-1][K-1];
            }
        }
        for (int J = JC1; J <= JB; J++) {
            for (int K = 3; K <= 10; K++) {
                int KK = IMAX[J-1] + 1;
                A[KK-1][J-1][K-1] = A[KK-1][J - 1-1][K-1];
            }
        }
        double DX2;
        double WK;
        for (int J = JC1; J <= JB; J++) {
            DX2 = X2[J-1] - X2[J-1 - 1];
            int KK = IMAX[J-1] + 1;
            WK = 1.0;
            WX = A[KK-1][J - 1-1][NF-1] / (DX2 * DX2);
            WY = A[KK-1][J - 1-1][NW-1] * ((R[J-1] - WK * DX2) * ((A[KK-1][J-1][NRO3-1] * (R[J-1] / 6.0 - WK * DX2 / 8.0)) + (A[KK-1][J - 1-1][NRO3-1] - A[KK-1][J-1][NRO3-1]) * (R[J-1] / 8.0 - WK * DX2 / 10.0)));
            WZ = R[J-1] * ((A[KK-1][J-1][NRO3-1] * (R[J-1] / 3.0 - 5.0 * WK * DX2 / 24.0)) + ((A[KK-1][J - 1-1][NRO3-1] - A[KK-1][J-1][NRO3-1]) * (5.0 * R[J-1] / 24.0 - 3.0 * WK * DX2 / 20.0)));
            A[KK-1][J-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (NITER <= NPRINT) {
                A[KK-1][J-1][NW-1] = 0.0;
            }
            if (A[KK-1][J-1][NW-1] >= 1.0E+14) {
                A[KK-1][J-1][NW-1] = 1.0E+14;
            }
            if (A[KK-1][J-1][NW-1] <= -1.0E+14) {
                A[KK-1][J-1][NW-1] = -1.0E+14;
            }
        }
        double DX1 = X1[IC2-1] - X1[IC2 - 1-1];
        for (int I = 2; I <= IC2; I++) {
            for (int K = 3; K <= 10; K++) {
                A[I-1][JA3-1][K-1] = A[I-1][JA3 - 1-1][K-1];
            }
        }
        DX2 = X2[JA3-1] - X2[JA3 - 1-1];
        WK = 1.0;
        int IC17 = IC2 - 1;
        for (int I = 2; I <= IC17; I++) {
            Z = A[I-1][JA3-1][NW-1];
            WX = (A[I-1][JA3 - 1-1][NF-1] - A[I-1][JA3-1][NF-1]) / (DX2 * DX2);
            WY = A[I-1][JA3 - 1-1][NW-1] * (R[JA3-1] - WK * DX2) * (A[I-1][JA3-1][NRO3-1] * (R[JA3-1] / 6.0 - WK * DX2 / 8.0) + (A[I-1][JA3 - 1-1][NRO3-1] - A[I-1][JA3-1][NRO3-1]) * (R[JA3-1] / 8.0 - WK * DX2 / 10.0));
            WZ = R[JA3-1] * (A[I-1][JA3-1][NRO3-1] * (R[JA3-1] / 3.0 - 5.0 * WK * DX2 / 24.0) + (A[I-1][JA3 - 1-1][NRO3-1] - A[I-1][JA3-1][NRO3-1]) * (5.0 * R[JA3-1] / 24.0 - 3.0 * WK * DX2 / 20.0));
            A[I-1][JA3-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (NITER <= NPRINT) {
                A[I-1][JA3-1][NW-1] = 0.0;
            }
            if (A[I-1][JA3-1][NW-1] > 1.0E+14) {
                A[I-1][JA3-1][NW-1] = 1.0E+14;
            }
            if (A[I-1][JA3-1][NW-1] < -1.0E+14) {
                A[I-1][JA3-1][NW-1] = -1.0E+14;
            }
        }
        int IC15 = IC2 + 1;
        int IC57 = IC1;
        for (int I = IC15; I <= IC57; I++) {
            for (int K = 3; K <= 10; K++) {
                A[I-1][JN-1][K-1] = A[I-1][JN - 1-1][K-1];
            }
        }
        DX2 = X2[JB1-1] - X2[JB1 - 1-1];
        WK = 1.0;
        for (int I = IC15; I <= IC57; I++) {
            Z = A[I-1][JB1-1][NW-1];
            WX = (A[I-1][JB1 - 1-1][NF-1] - A[I-1][JB1-1][NF-1]) / (DX2 * DX2);
            WY = A[I-1][JB1 - 1-1][NW-1] * (R[JB1-1] - WK * DX2) * (A[I-1][JB1-1][NRO3-1] * (R[JB1-1] / 6.0 - WK * DX2 / 8.0) + (A[I-1][JB1 - 1-1][NRO3-1] - A[I-1][JB1-1][NRO3-1]) * (R[JB1-1] / 8.0 - WK * DX2 / 10.0));
            WZ = R[JB1-1] * (A[I-1][JB1-1][NRO3-1] * (R[JB1-1] / 3.0 - 5.0 * WK * DX2 / 24.0) + (A[I-1][JB1 - 1-1][NRO3-1] - A[I-1][JB1-1][NRO3-1]) * (5.0 * R[JB1-1] / 24.0 - 3.0 * WK * DX2 / 20.0));
            A[I-1][JB1-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (NITER <= NPRINT) {
                A[I-1][JB1-1][NW-1] = 0.0;
            }
            if (A[I-1][JB1-1][NW-1] >= 1.0E+14) {
                A[I-1][JB1-1][NW-1] = 1.0E+14;
            }
            if (A[I-1][JB1-1][NW-1] <= -1.0E+14) {
                A[I-1][JB1-1][NW-1] = -1.0E+14;
            }
        }


        for (int I = IC; I <= INM; I++) {
            for (int K = 3; K <= 10; K++) {
                A[I-1][JC-1][K-1] = A[I-1][JC - 1-1][K-1];
            }
        }

        WK = 1.0;
        DX2 = X2[JC-1] - X2[JC -1- 1];
        for (int I = IC; I <= IN; I++) {
            Z = A[I-1][JC-1][NW-1];
            WX = A[I-1][JC - 1-1][NF-1] / (DX2 * DX2);
            WY = A[I-1][JC - 1-1][NW-1] * (R[JC-1] - WK * DX2) * (A[I-1][JC-1][NRO3-1] * (R[JC-1] / 6.0 - WK * DX2 / 8.0)
                    + (A[I-1][JC -1- 1][NRO3-1] - A[I-1][JC-1][NRO3-1]) * (R[JC-1] / 8.0 - WK * DX2 / 10.0));
            WZ = R[JC-1] * (A[I-1][JC-1][NRO3-1] * (R[JC-1] / 3.0 - 5.0 * WK * DX2 / 24.0)
                    + (A[I-1][JC -1- 1][NRO3-1] - A[I-1][JC-1][NRO3-1]) * (5.0 * R[JC-1] / 24.0 - 3.0 * WK * DX2 / 20.0));
            A[I-1][JC-1][NW-1] = -1.0 * (WX + WY) / WZ;

            if (NITER <= NPRINT) {
                A[I-1][JC-1][NW-1] = 0.0;
            }
            if (A[I-1][JC-1][NW-1] >= 1.0E+14) {
                A[I-1][JC-1][NW-1] = 1.0E+14;
            }
            if (A[I-1][JC-1][NW-1] <= -1.0E+14) {
                A[I-1][JC-1][NW-1] = -1.0E+14;
            }
        }

        for (int I = 2; I <= INM; I++) {
            for (int K = 3; K <= 10; K++) {
                A[I-1][1-1][K-1] = A[I-1][2-1][K-1];
                if (NCORD == 1) {
                    A[I-1][1-1][NW-1] = 0.0;
                }
                if (NCORD == 2) {
                    A[I-1][1-1][NW-1] = A[I-1][2-1][NW-1] + (X2[2-1] / (X2[3-1] - X2[2-1])) * (A[I-1][2-1][NW-1] - A[I-1][3-1][NW-1]);
                }
                if (A[I-1][1-1][NW-1] >= 1.0E+14) {
                    A[I-1][1-1][NW-1] = 1.0E+14;
                }
                if (A[I-1][1-1][NW-1] <= -1.0E+14) {
                    A[I-1][1-1][NW-1] = -1.0E+14;
                }
            }
        }
        for (int K = 1; K <= N3; K++) {
            for (int J = 1; J <= JC; J++) {
                A[IN-1][J-1][K-1] = A[IN -1- 1][J-1][K-1];
            }
        }
    }

    private static double sorce(int N1, int N2, int N3, double[][][] A, double SOURCE, int I, int J, int K) {
        switch (K) {
            case 1 -> SOURCE = UKS[I - 1][J - 1];
            case 2 -> SOURCE = A[I - 1][J - 1][NW - 1];
            case 3 -> SOURCE = 0;
            case 4 -> SOURCE = 0;
            case 5 -> SOURCE = 0;
            case 6 -> SOURCE = 0;
            default -> SOURCE = 0;
        }
        return SOURCE;
    }

    private static void convec(int N1, int N2, int N3, double[][][] A) {
        RL = 1.0;
        for (int J = 2; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) IL = IC2 + 1;
            if (J == JC) IM = IC - 1;
            for (int I = IL; I <= IM; I++) {
                double DV = R[J - 1] * (X1[I + 1 - 1] - X1[I - 1 - 1]) * (X2[J + 1 - 1] - X2[J - 1 - 1]);
                double G1PW = (A[I - 1][J + 1 - 1][NF - 1] - A[I - 1][J - 1 - 1][NF - 1] + A[I - 1 - 1][J + 1 - 1][NF - 1] - A[I - 1 - 1][J - 1 - 1][NF - 1]) / DV;
                double G1PE = (A[I - 1][J + 1 - 1][NF - 1] - A[I - 1][J - 1 - 1][NF - 1] + A[I + 1 - 1][J + 1 - 1][NF - 1] - A[I + 1 - 1][J - 1 - 1][NF - 1]) / DV;
                double G2PS = (A[I - 1 - 1][J - 1][NF - 1] - A[I - 1 + 1][J - 1][NF - 1] + A[I - 1 - 1][J - 1 - 1][NF - 1] - A[I + 1 - 1][J - 1 - 1][NF - 1]) / DV;
                double G2PN = (A[I - 1 - 1][J - 1][NF - 1] - A[I - 1 + 1][J - 1][NF - 1] + A[I - 1 - 1][J + 1 - 1][NF - 1] - A[I + 1 - 1][J + 1 - 1][NF - 1]) / DV;
                AE[I - 1][J - 1] = 0.5 * (Math.abs(G1PE) - G1PE);
                AW[I - 1][J - 1] = 0.5 * (Math.abs(G1PW) + G1PW);
                AN[I - 1][J - 1] = 0.5 * (Math.abs(G2PN) - G2PN);
                AS[I - 1][J - 1] = 0.5 * (Math.abs(G2PS) + G2PS);
                double BPP = A[I - 1][J - 1][NMU - 1];
                double BBE = (A[I + 1 - 1][J - 1][NMU - 1] + BPP) * BE[I - 1];
                double BBW = (A[I - 1 - 1][J - 1][NMU - 1] + BPP) * BW[I - 1];
                double BBN = (A[I - 1][J + 1 - 1][NMU - 1] + BPP) * BN[J - 1];
                double BBS = (A[I - 1][J - 1 - 1][NMU - 1] + BPP) * BS[J - 1];
                AE1[I - 1][J - 1] = AE[I - 1][J - 1] + BBE;
                AW1[I - 1][J - 1] = AW[I - 1][J - 1] + BBW;
                AN1[I - 1][J - 1] = AN[I - 1][J - 1] + BBN;
                AS1[I - 1][J - 1] = AS[I - 1][J - 1] + BBS;
            }
        }
    }

    private static void veldis(int N1, int N2, int N3, double[][][] A) {
        double H2;
        double RX2;
        double H1;
        double RX1;
        for (int J = 2; J <= JNM; J++) {
            H2 = (X2[J - 1] - X2[J - 1 - 1]) / (X2[J + 1 - 1] - X2[J - 1]);
            RX2 = R[J - 1] * (X2[J + 1 - 1] - X2[J - 1 - 1]);
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) IL = IC2 + 1;
            if (J == JC) IM = IC - 1;
            for (int I = IL; I <= IM; I++) {
                H1 = (X1[I - 1 - 1] - X1[I - 1]) / (X1[I + 1 - 1] - X1[I - 1]);
                RX1 = R[J - 1] * (X1[I + 1 - 1] - X1[I - 1 - 1]);
                A[I - 1][J - 1][NV1 - 1] = (A[I - 1][J + 1 - 1][NF - 1] - A[I - 1][J - 1][NF - 1]) * H2 + (A[I - 1][J - 1][NF - 1] - A[I - 1][J - 1 - 1][NF - 1]) / H2;
                A[I - 1][J - 1][NV1 - 1] = A[I - 1][J - 1][NV1 - 1] / RX2 / A[I - 1][J - 1][NRO3 - 1];
                A[I - 1][J - 1][NV2 - 1] = (A[I + 1 - 1][J - 1][NF - 1] - A[I - 1][J - 1][NF - 1]) * H1 + (A[I - 1][J - 1][NF - 1] - A[I - 1 - 1][J - 1][NF - 1]) / H1;
                A[I - 1][J - 1][NV2 - 1] = A[I - 1][J - 1][NV2 - 1] / RX1 / A[I - 1][J - 1][NRO3 - 1];
                if (NITER >= 85 && Math.abs(A[I - 1][J - 1][NV1 - 1]) >= 27000.0) NITER = NITER + 1000;
            }
        }
        double DELX2 = X2[1] - X2[0];
        for (int I = 2; I <= INM; I++) {
            A[I - 1][1 - 1][NV1 - 1] = (A[I - 1][2 - 1][NF - 1] - A[I - 1][1 - 1][NF - 1]) / DELX2 / A[I - 1][1 - 1][NRO3 - 1] / 0.5 / R[2 - 1];
        }
        for (int J = 2; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) IL = IC2 + 1;
            if (J == JC) IM = IC - 1;
            for (int I = IL; I <= IM; I++) {
                VS[I - 1][J - 1] = 0.5 * (A[I - 1][J - 1][NV1 - 1] * A[I - 1][J - 1][NV1 - 1] + A[I - 1][J - 1][NV2 - 1] * A[I - 1][J - 1][NV2 - 1]);
                FKS[I - 1][J - 1] = Math.pow(A[I - 1][J - 1][NMU2 - 1] / R[J - 1], 2) * A[I - 1][J - 1][NRO3 - 1];
            }
        }
        double DIF1;
        double DIF2;
        double DIF3;
        double DIF4;
        double DIF5;
        double S1;
        double S2;
        for (int J = 2; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JA3) IL = IC2 + 1;
            if (J == JC) IM = IC - 1;
            for (int I = IL; I <= IM; I++) {
                DIF1 = ((VS[I + 1 - 1][J - 1] - VS[I - 1][J - 1]) * (X1[I - 1] - X1[I - 1 - 1]) /
                        (X1[I + 1 - 1] - X1[I - 1]) + (VS[I - 1][J - 1] - VS[I - 1 - 1][J - 1]) * (X1[I + 1 - 1] - X1[I - 1]) /
                        (X1[I - 1] - X1[I - 1 - 1])) / (X1[I + 1 - 1] - X1[I - 1 - 1]);
                DIF2 = ((VS[I - 1][J + 1 - 1] - VS[I - 1][J - 1]) * (X2[J - 1] - X2[J - 1 - 1]) /
                        (X2[J + 1 - 1] - X2[J - 1]) + (VS[I - 1][J - 1] - VS[I - 1][J - 1 - 1]) * (X2[J + 1 - 1] - X2[J - 1]) /
                        (X2[J - 1] - X2[J - 1 - 1])) / (X2[J + 1 - 1] - X2[J - 1 - 1]);
                DIF3 = ((A[I + 1 - 1][J - 1][NRO3 - 1] - A[I - 1][J - 1][NRO3 - 1]) * (X1[I - 1] - X1[I - 1 - 1]) /
                        (X1[I + 1 - 1] - X1[I - 1]) + (A[I - 1][J - 1][NRO3 - 1] - A[I - 1 - 1][J - 1][NRO3 - 1]) * (X1[I + 1 - 1] - X1[I - 1]) /
                        (X1[I - 1] - X1[I - 1 - 1])) / (X1[I + 1 - 1] - X1[I - 1 - 1]);
                DIF4 = ((A[I - 1][J + 1 - 1][NRO3 - 1] - A[I - 1][J - 1][NRO3 - 1]) * (X2[J - 1] - X2[J - 1 - 1]) /
                        (X2[J + 1 - 1] - X2[J - 1]) + (A[I - 1][J - 1][NRO3 - 1] - A[I - 1][J - 1 - 1][NRO3 - 1]) * (X2[J + 1 - 1] - X2[J - 1]) /
                        (X2[J - 1] - X2[J - 1 - 1])) / (X2[J + 1 - 1] - X2[J - 1 - 1]);
                DIF5 = 0.0;
                DIF5 = ((FKS[I + 1 - 1][J - 1] - FKS[I - 1][J - 1]) * (X1[I - 1] - X1[I - 1 - 1]) /
                        (X1[I + 1 - 1] - X1[I - 1]) + (FKS[I - 1][J - 1] - FKS[I - 1 - 1][J - 1]) * (X1[I + 1 - 1] - X1[I - 1]) /
                        (X1[I - 1] - X1[I - 1 - 1])) / (X1[I + 1 - 1] - X1[I - 1 - 1]);
                S1 = DIF1 * DIF4;
                S2 = DIF2 * DIF3;
                UKS[I - 1][J - 1] = R[J - 1] * (S1 - S2) - DIF5;
            }
        }
    }

    private static void viscos(int N1, int N2, int N3, double[][][] A) {
        RL = 0.1;
        for (int K = 3; K <= 9; K++) { //to 100
            for (int J = 1; J <= JN; J++) { //to 100
                int IL = IMIN[J - 1];
                int IM = IMAX[J - 1] + 1;
                if (J <= JC) {
                    IM = IN;
                }
                if (J < JA1 && J > JA) {
                    IL = IMIN[J - 1] - 1;
                }
                if (J > JA2) {
                    IL = IMIN[J - 1] - 1;
                }
                for (int I = IL; I <= IM; I++) {
                    if (A[I - 1][J - 1][K - 1] < 0) {
                        A[I - 1][J - 1][K - 1] = 0;
                    }
                }
            }
        }
        double CM;
        for (int J = 1; J <= JN; J++) { //to 121
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J < JA1 && J > JA) {
                IL = IMIN[J - 1] - 1;
            }
            if (J < JZ1 && J > JZ) {
                IL = IMIN[J - 1] - 1;
            }
            if (J > JA2) {
                IL = IMIN[J - 1] - 1;
            }
            for (int I = IL; I <= IM; I++) {
                CM = A[I - 1][J - 1][NMFU1 - 1] + A[I - 1][J - 1][NMFU2 - 1] + A[I - 1][J - 1][NMOX1 - 1] + A[I - 1][J - 1][NMOX2 - 1] + A[I - 1][J - 1][NMPR1 - 1] + A[I - 1][J - 1][NMPR2 - 1] + A[I - 1][J - 1][NMDF - 1];
                if (CM <= 0) {
                    continue;
                }
                A[I - 1][J - 1][NMFU1 - 1] = A[I - 1][J - 1][NMFU1 - 1] / CM;
                A[I - 1][J - 1][NMFU2 - 1] = A[I - 1][J - 1][NMFU2 - 1] / CM;
                A[I - 1][J - 1][NMOX1 - 1] = A[I - 1][J - 1][NMOX1 - 1] / CM;
                A[I - 1][J - 1][NMOX2 - 1] = A[I - 1][J - 1][NMOX2 - 1] / CM;
                A[I - 1][J - 1][NMPR1 - 1] = A[I - 1][J - 1][NMPR1 - 1] / CM;
                A[I - 1][J - 1][NMPR2 - 1] = A[I - 1][J - 1][NMPR2 - 1] / CM;
                A[I - 1][J - 1][NMDF - 1] = A[I - 1][J - 1][NMDF - 1] / CM;
            }
        } //RK4 48
        double AZK;
        double ZOKD;
        double ZGD;
        double XXX;
        double AA = 0;
        double BA = 0;
        double CA = 0;
        double DA = 0;
        double EA = 0;
        double YYY;
        double Z;
        double CK;
        double CG;
        double TQ;
        double TZ;
        double DY1;
        double DYY;
        double ROQ;
        double ROZ;
        double ROX;
        double RK4;
        double VK1;
        double RZ1;
        double RK;
        double RK1;
        double RK2 = 0;
        for (int J = 1; J <= JN; J++) { //to 20
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IL;
            }
            if (J < JA1 && J > JA) {
                IL = IMIN[J - 1] - 1;
            }
            if (J < JZ1 && J > JZ) {
                IL = IMIN[J - 1] - 1;
            }
            if (J > JA2) {
                IL = IMIN[J - 1] - 1;
            }
            for (int I = IL; I <= IM; I++) {
                AZK = A[I - 1][J - 1][NMPR1 - 1] + A[I - 1][J - 1][NMOX2 - 1] + A[I - 1][J - 1][NMDF - 1] + A[I - 1][J - 1][NMPR2 - 1] + A[I - 1][J - 1][NMFU2 - 1];
                ZOKD = AZK * 0.7966 + A[I - 1][J - 1][NMOX1 - 1];
                ZGD = AZK * 0.2034 + A[I - 1][J - 1][NMFU1 - 1];
                if (ZGD < 0.005) {
                    ZGD = 0.005;
                }
                XXX = ZOKD / ZGD / 3.9166;
                if (XXX <= 0.8) {
                    AA = 0.02323446;
                    BA = 5968.083;
                    CA = 9185.521;
                    DA = -24199.95;
                }
                if (XXX > 0.8) {
                    if (A[I - 1][J - 1][NK - 1] <= 2) {
                        AA = -0.3167559;
                        BA = 10790.7;
                        CA = -11990.44;
                        DA = 5632.969;
                        EA = -969.617;
                        //goto 2
                    }
                    if (XXX > 2) {
                        AA = -1.954687;
                        BA = 4416.632;
                        CA = -2128.663;
                        DA = 412.6477;
                        EA = -29.03693;
                        //goto 5
                    }
                }

                // 5
                YYY = AA + BA * XXX + CA * XXX * XXX + DA * XXX * XXX * XXX + EA * XXX * XXX * XXX * XXX;
                if (YYY < 293) {
                    YYY = 293;
                }
                Z = A[I - 1][J - 1][NT - 1];
                CK = 0.3868;
                CG = 3.332;
                TQ = 3.332 * A[I - 1][J - 1][NMFU1 - 1] + 14.743 * A[I - 1][J - 1][NMDF - 1] + 1.081 * A[I - 1][J - 1][NMPR1 - 1] + 0.3868 * A[I - 1][J - 1][NMOX1 - 1] + 2.059 * A[I - 1][J - 1][NMPR2 - 1] + 1.105 * A[I - 1][J - 1][NMOX2 - 1] + 1.753 * A[I - 1][J - 1][NMFU2 - 1];
                TZ = A[I - 1][J - 1][NMU1 - 1] + 4681.06 * A[I - 1][J - 1][NMFU1 - 1] - A[I - 1][J - 1][NMFU2 - 1] * 2299.57 + A[I - 1][J - 1][NMOX2 - 1] * 3949.64 + 8950. * A[I - 1][J - 1][NMPR1 - 1] + 13445.55 * A[I - 1][J - 1][NMPR2 - 1];
                if (TZ < 0) {
                    TZ = 0;
                }
                if (TQ > 1.E-10) {
                    A[I - 1][J - 1][NT - 1] = TZ / TQ + 293; // RK4 97
                    if (XXX > 0.8) {
                        A[I - 1][J - 1][NT - 1] = A[I - 1][J - 1][NT - 1] * 1.1;
                    }
                    if (XXX >= 0.2 && XXX <= 0.8) {
                        A[I - 1][J - 1][NT - 1] = A[I - 1][J - 1][NT - 1] * 1.2;
                    }
                    //goto 50
                }
                //50
                if (A[I - 1][J - 1][NT - 1] < 293) {
                    A[I - 1][J - 1][NT - 1] = 293;
                }
                if (A[I - 1][J - 1][NT - 1] < 4500) {
                    A[I - 1][J - 1][NT - 1] = 4500;
                }
                DY1 = A[I - 1][J - 1][NT - 1] - YYY;
                DYY = (A[I - 1][J - 1][NT - 1] - YYY) / A[I - 1][J - 1][NT - 1];
                Z = A[I][J][NRO3];
                ROQ = A[I - 1][J - 1][NT - 1] * (A[I - 1][J - 1][NMFU1 - 1] / 16. + A[I - 1][J - 1][NMOX1 - 1] / 32. + A[I - 1][J - 1][NMDF - 1] / 2. + A[I - 1][J - 1][NMOX2 - 1] / 28. + A[I - 1][J - 1][NMPR1 - 1] / 44. + A[I - 1][J - 1][NMPR2 - 1] / 18. + A[I - 1][J - 1][NMFU2 - 1] / 17.);
                if (ROQ != 0) {
                    ROZ = A[I - 1][J - 1][NMOX1 - 1] + A[I - 1][J - 1][NMOX2 - 1] + A[I - 1][J - 1][NMDF - 1] + A[I - 1][J - 1][NMFU1 - 1] + A[I - 1][J - 1][NMPR1 - 1] + A[I - 1][J - 1][NMPR2 - 1] + A[I - 1][J - 1][NMFU2 - 1];
                    ROX = PKS[I - 1] / 8314.4;
                    A[I - 1][J - 1][NRO3 - 1] = ROX * ROZ / ROQ;
                    A[I - 1][J - 1][NRO3 - 1] = Z + (A[I - 1][J - 1][NRO3 - 1] - Z) * RL;
                    //goto 60
                }
                // 60
                if (A[I - 1][J - 1][NRO3 - 1] < 0.00001) {
                    A[I - 1][J - 1][NRO3 - 1] = 0.00001;
                }
                // 20 RK4 118
            }// 20
        }// 20
        RK4 = 0;
        for (int J = 2; J <= JA3; J++) {
            RK4 = RK4 + 3.1416 * (A[IC2 - 1][J - 1][NF - 1] - A[IC2 - 1][J - 1 - 1][NF - 1]) * (Math.pow(A[IC2 - 1][J - 1][NV1 - 1], 2) + Math.pow(A[IC2 - 1][J - 1][NMU2 - 1], 2));
        }
        for (int J = 1; J <= JN; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J == JC) {
                IM = IN;
            }
            for (int I = IL; I <= IM; I++) {
                VK1 = Math.pow((VINP * VINP + Math.pow((A[1 - 1][JA - 1][NMU2 - 1] / R[JA - 1 - 1]), 2)), 0.5);
                RZ1 = Math.pow((2. * R[JN - 1]), 0.666) / Math.pow((X1[IC1 - 1]), 0.333);
                if (I <= IC2) {
                    RZ1 = Math.pow((2 * R[JA3 - 1]), 0.666) / Math.pow((X1[IC1 - 1]), 0.333);
                }
                RK = 3.1416 * (R[JA2 - 1] * R[JA2 - 1] - R[JA1 - 1] * R[JA1 - 1]) * ROS * Math.pow(VINS, 3);
                RK1 = 3.1416 * (R[JA - 1] - R[JZ1 - 1]) * (R[JA - 1] + R[JZ1 - 1]) * ROP * VINP * VINP * VINP + 3.1416 * R[JZ - 1] * R[JZ - 1] * A[1 - 1][2 - 1][NRO3 - 1] * Math.pow(VINS1, 3);
                if (I > IC2) {
                    RK2 = RK + RK4;
                }
                if (I <= IC2) {
                    RK2 = RK1;
                }
                A[I - 1][J - 1][NMU - 1] = 0.012 * RZ1 * Math.pow(A[I - 1][J - 1][NRO3 - 1], 0.666) * Math.pow(RK2, 0.333);
                A[I - 1][J - 1][NMU - 1] = A[I - 1][J - 1][NMU - 1] * 2;
            }
        }
    }

    private static void grid(int N1, int N2, int N3, double[][][] A) {
        int NIG = 0;
        NIG = NIG + 1;
        X2[JZ - 1] = DCMB / 2.0;
        X2[JZ1 - 1] = DOGOR / 2.0;
        X2[JA - 1] = DELGOR + X2[JZ1 - 1];
        X2[JA3 - 1] = DKAR / 2.0;
        X2[JA1 - 1] = DOOK / 2.0;
        X2[JA2 - 1] = DELOK + X2[JA1 - 1];
        X2[JN - 1] = DKS / 2.0;
        X1[IC2 - 1] = DELKAR;
        X1[IC1 - 1] = X1[IC2 - 1] + ALKS;
        X1[IC - 1] = ALKD;
        for (int J = 2; J <= JZ; J++) {
            X2[J - 1] = X2[J - 1 - 1] + X2[JZ - 1] / (JZ - 1);
        }
        int JZZ = JZ + 1;
        for (int J = JZZ; J <= JZ1; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JZ1 - 1] - X2[JZ - 1]) / (JZ1 - JZ);
        }
        JZZ = JZ1 + 1;
        for (int J = JZZ; J <= JA; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JA - 1] - X2[JZ1 - 1]) / (JA - JZ1);
        }
        JZZ = JA + 1;
        for (int J = JZZ; J <= JA3; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JA3 - 1] - X2[JA - 1]) / (JA3 - JA);
        }
        JZZ = JA3 + 1;
        for (int J = JZZ; J <= JA1; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JA1 - 1] - X2[JA3 - 1]) / (JA1 - JA3);
        }
        JZZ = JA1 + 1;
        for (int J = JZZ; J <= JA2; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JA2 - 1] - X2[JA1 - 1]) / (JA2 - JA1);
        }
        JZZ = JA2 + 1;
        for (int J = JZZ; J <= JN; J++) {
            X2[J - 1] = X2[J - 1 - 1] + (X2[JN - 1] - X2[JA2 - 1]) / (JN - JA2);
        }
        int KJ;
        int JJ = 0;
        for (int J = 2; J <= JN; J++) {
            double XX2 = DKR / 2.0 - X2[J];
            JJ = J;
            if (XX2 <= 0.0) break;
        }
        KJ = JJ - JC;
        JC = JJ;
        IC = IC - KJ;
        int JB1 = JN;
        int JB = JB1 - 1;
        int JC1 = JC + 1;
        X2[JC - 1] = DKR / 2.0;
        int IC1 = IC - (JN - JC);
        for (int I = 2; I <= IC2; I++) {
            X1[I - 1] = X1[I - 1 - 1] + X1[IC2 - 1] / (IC2 - 1);
        }
        int ICC = IC2 + 1;
        X1[IC1 - 1] = X1[IC2 - 1] + ALKS;
        for (int I = ICC; I <= IC1; I++) {
            X1[I - 1] = X1[I - 1 - 1] + (X1[IC1 - 1] - X1[IC2 - 1]) / (IC1 - IC2);
        }
        int ICLM1 = IC1 + 1;
        for (int I = ICLM1; I <= IC; I++) {
            int LK = JN - I + IC1;
            X1[I - 1] = X1[IC1 - 1] + (X2[JN - 1] - X2[LK - 1]) / TGY;
        }
        ICC = IC + 1;
        for (int I = ICC; I <= IN; I++) {
            double XZ1 = 1.0;
            if (NVAR != 1) XZ1 = 1.0E-3;
            X1[I - 1] = X1[I - 1 - 1] + XZ1 / (IN - IC);
        }
        if (NVAR == 1) {
            for (int J = 1; J <= JN; J++) {
                X2[J - 1] = X2[J - 1] * 0.001;
            }
            for (int I = 1; I <= IN; I++) {
                X1[I - 1] = X1[I - 1] * 0.001;
            }
        }
        if (NCORD == 14 || NCORD == 1) {
            for (int J = 1; J <= JN; J++) {
                R[J - 1] = 1.0;
            }
        } else if (NCORD == 15 || NCORD == 2) {
            for (int J = 1; J <= JN; J++) {
                R[J - 1] = X2[J - 1];
            }
        }
        for (int I = 2; I <= INM; I++) {
            double DX1 = 1.0 / (X1[I + 1 - 1] - X1[I - 1 - 1]);
            BW[I - 1] = DX1 / (X1[I - 1] - X1[I - 1 - 1]);
            BE[I - 1] = DX1 / (X1[I + 1 - 1] - X1[I - 1]);
        }
        for (int J = 2; J <= JNM; J++) {
            double DX2 = 0.5 / (X2[J + 1 - 1] - X2[J - 1 - 1]);
            BS[J - 1] = (1.0 + R[J - 1 - 1] / R[J - 1]) / (X2[J - 1] - X2[J - 1 - 1]) * DX2;
            BN[J - 1] = (1.0 + R[J + 1 - 1] / R[J - 1]) / (X2[J - 1 + 1] - X2[J - 1]) * DX2;
        }
        for (int J = 1; J <= JN; J++) {
            for (int I = 1; I <= IN; I++) {
                AE[I - 1][J - 1] = 0.0;
                AW[I - 1][J - 1] = 0.0;
                AN[I - 1][J - 1] = 0.0;
                AS[I - 1][J - 1] = 0.0;
                AE1[I - 1][J - 1] = 0.0;
                AW1[I - 1][J - 1] = 0.0;
                AN1[I - 1][J - 1] = 0.0;
                AS1[I - 1][J - 1] = 0.0;
                UKS[I - 1][J - 1] = 0.0;
                VS[I - 1][J - 1] = 0.0;
                FKS[I - 1][J - 1] = 0.0;
                TKS[I - 1][J - 1] = 0.0;
            }
        }
        for (int K = 1; K <= N3; K++) {
            for (int J = 1; J <= JN; J++) {
                for (int I = 1; I <= IN; I++) {
                    A[I - 1][J - 1][K - 1] = 0.0;
                }
            }
        }
    }

    public static void write7(String str) {
        System.out.println(str);
    }

    public static void write9(String str) {
        System.out.println(str);
    }

    public static void init(int N1, int N2, int N3, double[][][] A) {
        for (int J = 1; J <= JA3; J++) {
            IMIN[J - 1] = 2;
        }
        int JAA = JA3 + 1;
        for (int J = JAA; J <= JN; J++) {
            IMIN[J - 1] = IN - 1;
        }
        for (int J = JC1; J <= JN; J++) {
            IMAX[J - 1] = IC - 1 + JC - J;
        }
        for (int J = JA1; J <= JA2; J++) {
            A[IC2 - 1][J - 1][NV1 - 1] = VINS;
            A[IC2 - 1][J - 1][NMFU1 - 1] = 0;
            A[IC2 - 1][J - 1][NRO3 - 1] = ROS;
            A[IC2 - 1][J - 1][NMOX1 - 1] = ZMO2N;
            A[IC2 - 1][J - 1][NT - 1] = 298;
            A[IC2 - 1][J - 1][NMU1 - 1] = 0;
            if (NCORD == 2) {
                A[IC2 - 1][J - 1][NF - 1] = ROS * VINS * Math.pow(R[J - 1] - R[JA2 - 1], 0.5);
            }
            if (NCORD == 1) {
                A[IC2 - 1][J - 1][NF - 1] = A[IC2 - 1][JA2 - 1][NF - 1] - ROS * VINS * (X2[JA2 - 1] - X2[J - 1]);
            }
        }
        A[IC2 - 1][JA2 - 1][NV1 - 1] = 0;
        for (int J = JA3; J <= JA1; J++) {
            A[IC2 - 1][J - 1][NF - 1] = A[IC2 - 1][JA1 - 1][NF - 1];
        } // RK1 91
        for (int I = 1; I <= IC2; I++) {
            A[I - 1][JA3 - 1][NF - 1] = A[IC2 - 1][JA3 - 1][NF - 1];
        }
        A[IC2 - 1][JA1 - 1][NV1 - 1] = 0;
        for (int J = JA; J <= JA3; J++) {
            A[1 - 1][J - 1][NF - 1] = A[1 - 1][JA3 - 1][NF - 1];
        }
        for (int J = JZ1; J <= JA; J++) {
            A[1 - 1][J - 1][NV1 - 1] = VINP;
            A[1 - 1][J - 1][NMU2 - 1] = VINP * Math.pow(R[J - 1], 2) / R[JA - 1] * 1 / 5;
            A[1 - 1][J - 1][NMOX1 - 1] = 0;
            A[1 - 1][J - 1][NMFU1 - 1] = ZMF2N;
            A[1 - 1][J - 1][NRO3 - 1] = ROP;
            A[1 - 1][J - 1][NT - 1] = 298;
            A[1 - 1][J - 1][NMFU1 - 1] = H;
            if (NCORD == 2) {
                A[1 - 1][J - 1][NF - 1] = A[1 - 1][JA - 1][NF - 1] + ROP * VINP * (R[J - 1] - R[JA - 1]) * 0.5 * (R[J - 1] + R[JA - 1]);
            }
            if (NCORD == 1) {
                A[1 - 1][J - 1][NMU1 - 1] = A[1 - 1][JA - 1][NF - 1] + ROP * VINP * (X2[J - 1] - X2[JA - 1]);
            }
        } //61
        A[1 - 1][JA - 1][NV1 - 1] = 0;
        A[1 - 1][JZ1 - 1][NV1 - 1] = 0;
        for (int J = JZ; J <= JZ1; J++) {
            A[1 - 1][J - 1][NT - 1] = A[1 - 1][JZ1 - 1][NF - 1];
        }
        for (int J = 1; J <= JZ; J++) {
            A[1 - 1][J - 1][NT - 1] = 298;
            A[1 - 1][J - 1][NV1 - 1] = VINS1;
            A[1 - 1][J - 1][NMOX1 - 1] = 1;
            A[1 - 1][J - 1][NRO3 - 1] = 1. / (A[1 - 1][J - 1][NMFU1 - 1] / 16. + A[1 - 1][J - 1][NMOX1 - 1] / 32. + A[1 - 1][J - 1][NMPR1 - 1] / 44. + A[1 - 1][J - 1][NMOX2 - 1] / 28. + A[1 - 1][J - 1][NMPR2 - 1] / 18. + A[1 - 1][J - 1][NMFU2 - 1] / 17. + A[1 - 1][J - 1][NMDF - 1] / 2);
            A[1 - 1][J - 1][NRO3 - 1] = P * A[1 - 1][J - 1][NRO3 - 1] / 8314.4 / A[1 - 1][J - 1][NT - 1];
            if (NCORD == 2) {
                A[1 - 1][J - 1][NF - 1] = A[1 - 1][JZ - 1][NF - 1] + VINS1 * A[1 - 1][J - 1][NRO3 - 1] * (R[J - 1] - R[JZ - 1]) * 0.5 * (R[J - 1] + R[JZ - 1]);
            }
        }
        for (int I = 1; I <= IN; I++) {
            A[I - 1][1 - 1][NF - 1] = A[1 - 1][1 - 1][NF - 1];
        }
        A[1 - 1][JZ - 1][NV1 - 1] = 0;
        for (int J = JZ1; J <= JA; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J == JC1 || J == JB) {
                IM = IMAX[J - 1] + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMFU1 - 1] = 1;
                A[I - 1][J - 1][NMOX1 - 1] = 0;
                A[I - 1][J - 1][NRO3 - 1] = ROP;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NMU1 - 1] = H;
            }
        } //RK1 145
        int JA17 = JA + 1;
        int JA18 = JA1 - 1;
        for (int J = JA17; J <= JA18; J++) {
            int IL = IMIN[J - 1] - 1;
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J == JC1 || J == JB) {
                IM = IMAX[J - 1] + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = 1;
                A[I - 1][J - 1][NMFU2 - 1] = 0;
                A[I - 1][J - 1][NRO3 - 1] = ROS;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NNFU1 - 1] = 0;
                A[I - 1][J - 1][NMU1 - 1] = H / 100;
            }
        }
        JA17 = JA1 + 1;
        for (int J = JA1; J <= JN; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J == JC1 || J == JB) {
                IM = IMAX[J - 1] + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMFU1 - 1] = 0;
                A[I - 1][J - 1][NT - 1] = 298;
                A[I - 1][J - 1][NRO3 - 1] = ROS;
                A[I - 1][J - 1][NMFU2 - 1] = 0;
                A[I - 1][J - 1][NMOX1 - 1] = 1;
                A[I - 1][J - 1][NMU1 - 1] = H / 100;
            }
        } // RK1 174
        for (int J = 1; J <= JZ; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J == JB || J == JC1) {
                IM = IMAX[J - 1] + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = A[1 - 1][J - 1][NMOX1 - 1];
                A[I - 1][J - 1][NRO3 - 1] = A[1 - 1][J - 1][NRO3 - 1];
                A[I - 1][J - 1][NMU1 - 1] = A[1 - 1][J - 1][NMU1 - 1];
                A[I - 1][J - 1][NMPR1 - 1] = A[1 - 1][J - 1][NMPR1 - 1];
                A[I - 1][J - 1][NMPR2 - 1] = A[1 - 1][J - 1][NMPR2 - 1];
                A[I - 1][J - 1][NMFU1 - 1] = A[1 - 1][J - 1][NMFU1 - 1];
                A[I - 1][J - 1][NT - 1] = ROM1 + 10;
            }
        }
        int JZ17 = JZ1 - 1;
        int JZ18 = JZ + 1;
        for (int J = JZ18; J <= JZ17; J++) {
            int IL = IMIN[J - 1] - 1;
            int IM = IMAX[J - 1] + 1;
            if (J <= JC) {
                IM = IN;
            }
            if (J == JB || J == JC1) {
                IM = IMAX[J - 1] + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NMOX1 - 1] = 0;
                A[I - 1][J - 1][NMFU1 - 1] = 1;
                A[I - 1][J - 1][NRO3 - 1] = ROP;
                A[I - 1][J - 1][NMU1 - 1] = H;
                A[I - 1][J - 1][NT - 1] = 298;
            }
        } // RK 1 202
        int JAB = JA3 - 1;
        for (int J = 2; J <= JAB; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JC) {
                IM = IC - 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NF - 1] = A[1 - 1][J - 1][NF - 1];
            }
        }
        for (int J = JA3; J <= JNM; J++) {
            int IL = IMIN[J - 1];
            int IM = IMAX[J - 1];
            if (J == JC) {
                IM = IC - 1;
            }
            if (J == JA3) {
                IL = IC2 + 1;
            }
            for (int I = IL; I <= IM; I++) {
                A[I - 1][J - 1][NF - 1] = A[IC2 - 1][J - 1][NF - 1];
            }
        }
    }

    public static void trueInit() {
        DKR = 8.0;
        DCMB = 2.0;
        DOGOR = 8.0;
        DELGOR = 1.0;
        DKAR = 12.0;
        DOOK = 20.0;
        DELOK = 0.4;
        DKS = 24.0;
        DELKAR = 10.0;
        ALKS = 15.0;
        ALKD = 42.239;
        TGY = 0.573;
        ALFA = 0.6;

        NW = 1;
        NF = 2;
        NMFU1 = 3;
        NMOX1 = 4;
        NMPR1 = 5;
        NMPR2 = 6;
        NMDF = 7;
        NMOX2 = 8;
        NMFU2 = 9;
        NMU1 = 10;
        NMU2 = 16;
        NRO3 = 12;
        NT = 13;
        NV1 = 14;
        NV2 = 15;
        NMU = 11;
        NK = 17;
        IE = 10;
        IV = 16;

        NMAX = 1000;
        NPRINT = 200;
        IP = 1;
        CC = 0.0035;
        NVARM = 5;

        VINP = 9.329155;
        VINS = 26.03549;
        ZMF2N = 1.00;
        ZMO2N = 1.00;
        VINT = 0.1;
        VINS1 = 2.E-3;
        ROP = 6.458;
        ROS = 12.916;

        RO1 = 6.458;
        ROM1 = 815.0;
        ROT = 0.1;

        NITER = 0;
        LN = 20;
        AFI = 0.34;
        AKT = 1.2;

        D = 5.E-3;
        P = 10.E+5;
        H = -4681.06;

        NCORD = 2;
        JA = 9;
        JA1 = 15;
        JC1 = 13;
        JB = 21;
        IC1 = 31;
        JA2 = 18;
        JB1 = 22;
        IC = 45;
        IN = 46;
        JN = 22;
        IC2 = 10;
        JZ = 4;
        JZ1 = 6;
        JA3 = 12;
        JC = 6;

        ROS1 = 6.458;
        ROP1 = 12.916;
        VINP1 = 0.1;
        VINT1 = 0.1;

        double[] X1Data = {0.0, 0.6, 1.2, 1.8, 2.4, 3.0, 4.4, 5.8, 7.2, 8.6, 10.0, 11.4, 12.8, 14.2, 15.6, 17.0, 18.4, 19.8, 21.2, 22.6, 24.0, 25.4, 26.8, 28.2, 29.6, 31.0, 32.4, 33.8, 35.2, 36.6, 38.0, 39.4, 40.8, 42.2, 43.6, 45.0, 46.4, 47.8, 49.2, 50.6, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0};
        double[] X2Data = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.6, 7.1, 8.0, 8.75, 9.5, 10.25, 11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0};

        System.arraycopy(X1Data, 0, X1, 0, X1.length);
        System.arraycopy(X2Data, 0, X2, 0, X2.length);

        System.arraycopy(new double[12], 0, RSDU, 0, RSDU.length);
        System.arraycopy(new int[25], 0, LMIN, 0, LMIN.length);
    }
}
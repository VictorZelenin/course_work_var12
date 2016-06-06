import mpi.MPI;

import java.util.Arrays;

public class Main {

    /**
     * quantity of nodes on x axis and t axis
     */
    private static final int QUANTITY_OF_X = 100;
    private static final int QUANTITY_OF_T = 32;

    /**
     * step on x-axis and t-axis
     */
    private static final double H = 1.0 / QUANTITY_OF_X;
    private static final double TAU = 1.0 / QUANTITY_OF_T;
    private static final double SIGMA = TAU / (H * H);

    /**
     * constants in function
     */
    private static final double a = 70;
    private static final double b = 1;
    private static final double C = 1.0;
    private static final double LAMBDA = 3;

    private static final double[] X_VECTOR;
    private static final double[] T_VECTOR;

    static {
        X_VECTOR = getXVector();
        T_VECTOR = getTVector();
    }

    private static double[] getXVector() {
        double[] xVector = new double[QUANTITY_OF_X + 1];

        for (int i = 0; i < xVector.length; i++) {
            xVector[i] = i * H;
        }

        return xVector;
    }

    private static double[] getTVector() {
        double[] tVector = new double[QUANTITY_OF_T + 1];

        for (int k = 0; k < tVector.length; k++) {
            tVector[k] = k * TAU;
        }

        return tVector;
    }

    private static double getExactSolution(double x, double t) {
        return Math.pow(C * Math.exp(-3.0 * LAMBDA / a * (x + LAMBDA * t) - b / (4.0 * LAMBDA)), -1.0 / 3.0);
    }
    
    // початкові умови
    private static double[] getStartVectorOmega() {
        double[] vector = new double[X_VECTOR.length];

        for (int i = 0; i < vector.length; i++) {
            vector[i] = getExactSolution(X_VECTOR[i], 0);
        }

        return vector;
    }

    // похідна по омега і - 1
    private static double getFirstMatrixCoefficient(double[] omega, int i) {
        return a * SIGMA - (b * TAU) / (2.0 * H) * Math.pow(omega[i], 3.0);
    }

    // похідна по омега і
    private static double getSecondMatrixCoefficient(double[] omega, int i) {
        return -2 * a * SIGMA + 3 * (b * TAU) / (2 * H) * Math.pow(omega[i], 2) * (omega[i + 1]
                - omega[i - 1]) - 1;
    }

    // похідна по омега і + 1
    private static double getThirdMatrixCoefficient(double[] omega, int i) {
        return a * SIGMA + (b * TAU) / (2 * H) * Math.pow(omega[i], 3);
    }

    // f(i)
    private static double getNonBoundaryFunctionValue(double[] oldOmega, double[] omega, int i) {
        return oldOmega[i] - omega[i] + a * SIGMA * (omega[i - 1] - 2 * omega[i] + omega[i + 1])
                + (b * TAU) / (2 * H) * Math.pow(omega[i], 3) * (omega[i + 1] - omega[i - 1]) / (2 * H);

    }

    // з лівого краю
    private static double getLeftBoundaryFunctionValue(double[] omega, int k) {
        return omega[0] - Math.pow(C * Math.exp(-3.0 * LAMBDA * LAMBDA / a * T_VECTOR[k]) - b / (4 * LAMBDA),
                -1.0 / 3.0);
    }

    // з правого краю
    public static double getRightBoundaryFunctionValue(double[] omega, int k) {
        return omega[QUANTITY_OF_X] - Math.pow(C * Math.exp(-3 * LAMBDA / a * (1 + LAMBDA * T_VECTOR[k])),
                -1.0 / 3);
    }

    // 3х діагональна матриця
    private static double[][] getCoefficientMatrix(double[] vectorOmega) {
        double[][] matrix = new double[QUANTITY_OF_X + 1][QUANTITY_OF_X + 1];

        // похідні з крайових умов
        matrix[0][0] = 1;
        matrix[QUANTITY_OF_X][QUANTITY_OF_X] = 1;

        for (int i = 1; i < QUANTITY_OF_X; i++) {
            matrix[i][i - 1] = getFirstMatrixCoefficient(vectorOmega, i);
            matrix[i][i] = getSecondMatrixCoefficient(vectorOmega, i);
            matrix[i][i + 1] = getThirdMatrixCoefficient(vectorOmega, i);
        }

        return matrix;
    }

    // вектор -F
    private static double[] getFunctionVector(double[] vectorOmegaOld, double[] vectorOmega, int k) {
        double[] vector = new double[QUANTITY_OF_X + 1];

        vector[0] = -getLeftBoundaryFunctionValue(vectorOmega, k);
        vector[QUANTITY_OF_X] = -getRightBoundaryFunctionValue(vectorOmega, k);


        for (int i = 1; i < QUANTITY_OF_X; i++) {
            vector[i] = -getNonBoundaryFunctionValue(vectorOmegaOld, vectorOmega, i);
        }

        return vector;
    }

    // w + delta w
    private static double[] getUpdatedVectorOmega(double[] vectorOmega, double[] vectorDeltaOmega) {
        double[] vector = new double[vectorOmega.length];

        for (int i = 0; i < vector.length; i++) {
            vector[i] = vectorOmega[i] + vectorDeltaOmega[i];
        }

        return vector;
    }

    private static double[][] calculateVectorOmega() {
        int processRank = MPI.COMM_WORLD.Rank();

        double[][] vectorOmegaInDifferentTimeLayers = new double[QUANTITY_OF_T + 1][QUANTITY_OF_X + 1];


        double[] coefficientArray = new double[(QUANTITY_OF_X + 1) * (QUANTITY_OF_X + 1)];
        double[] resultVector = new double[QUANTITY_OF_X + 1];

        if (processRank == 0) {
            vectorOmegaInDifferentTimeLayers[0] = getStartVectorOmega();
        }

        for (int k = 1; k < QUANTITY_OF_T + 1; k++) {
            if (processRank == 0) {
                double[][] coefficientMatrix = getCoefficientMatrix(vectorOmegaInDifferentTimeLayers[k - 1]);
                resultVector = getFunctionVector(vectorOmegaInDifferentTimeLayers[k - 1],
                        vectorOmegaInDifferentTimeLayers[k - 1], k);

                coefficientArray = Util.encode(coefficientMatrix);
            }

            MPI.COMM_WORLD.Bcast(coefficientArray, 0, coefficientArray.length, MPI.DOUBLE, 0);
            MPI.COMM_WORLD.Bcast(resultVector, 0, resultVector.length, MPI.DOUBLE, 0);

            PTMA.solveSystemOfEquations(coefficientArray, resultVector);

            if (processRank == 0) {
                double[] delta = new double[QUANTITY_OF_X + 1];
                MPI.COMM_WORLD.Recv(delta, 0, delta.length, MPI.DOUBLE, 1, 5);

                vectorOmegaInDifferentTimeLayers[k] =
                        getUpdatedVectorOmega(vectorOmegaInDifferentTimeLayers[k - 1], delta);

                System.out.print("on time [" + T_VECTOR[k] + "] х vector -> ");
                System.out.println(Arrays.toString(vectorOmegaInDifferentTimeLayers[k]));
            }
        }
        return vectorOmegaInDifferentTimeLayers;
    }

    private static double[] getRealVectorOmega(int k) {
        double[] vector = new double[X_VECTOR.length];

        for (int i = 0; i < vector.length; i++) {
            vector[i] = getExactSolution(X_VECTOR[i], T_VECTOR[k]);
        }

        return vector;
    }

    private static void calculateErrors(double[][] calculatedValues, double[][] realValues) {
        for (int i = 0; i < QUANTITY_OF_T; i++) {
            for (int j = 0; j <= QUANTITY_OF_X; j++) {
                double d = 100 * Math.abs(calculatedValues[i][j] - realValues[i][j]) / realValues[i][j];
                System.out.print(d + "% ");
            }

            System.out.println();
        }
    }

    public static void main(String[] args) {
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.Rank();

        double[][] calculated = calculateVectorOmega();

        if (rank == 0) {
            double[][] real = getRealOmegas();
            calculateErrors(calculated, real);
            Util.writeResultInFile("omegas.dat", calculated, X_VECTOR);
            Util.printVectorInFile("x-s.dat", X_VECTOR);
            Util.printVectorInFile("t-s.dat", T_VECTOR);
        }

        MPI.Finalize();
    }

    public static double[][] getRealOmegas() {
        double[][] matrix = new double[QUANTITY_OF_T + 1][QUANTITY_OF_X + 1];

        for (int i = 0; i < QUANTITY_OF_T + 1; i++) {
            for (int j = 0; j < QUANTITY_OF_X + 1; j++) {
                matrix[i][j] = getExactSolution(X_VECTOR[j], T_VECTOR[i]);
            }
        }

        return matrix;
    }


}
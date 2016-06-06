import mpi.MPI;

public class PTMA {

    public static void solveSystemOfEquations(double[] array, double[] b) {
        double[][] A = Util.decode(array);

        int n = b.length;
        int m = b.length / 2;

        int processRank = MPI.COMM_WORLD.Rank();

        double[] alpha = new double[m];
        double[] beta = new double[m];

        double[] epsilon = new double[n];
        double[] eta = new double[n];

        double[] x = new double[n];

        if (processRank == 0) {
            alpha = calculateAlpha(A, m);
            beta = calculateBeta(A, b, m, alpha);

            MPI.COMM_WORLD.Send(alpha, 0, alpha.length, MPI.DOUBLE, 1, 0);
            MPI.COMM_WORLD.Send(beta, 0, beta.length, MPI.DOUBLE, 1, 1);

            MPI.COMM_WORLD.Recv(epsilon, 0, epsilon.length, MPI.DOUBLE, 1, 2);
            MPI.COMM_WORLD.Recv(eta, 0, eta.length, MPI.DOUBLE, 1, 3);
        } else {
            epsilon = calculateEpsilon(A, n, m);
            eta = calculateEta(A, b, epsilon, n, m);

            MPI.COMM_WORLD.Recv(alpha, 0, alpha.length, MPI.DOUBLE, 0, 0);
            MPI.COMM_WORLD.Recv(beta, 0, beta.length, MPI.DOUBLE, 0, 1);

            MPI.COMM_WORLD.Send(epsilon, 0, epsilon.length, MPI.DOUBLE, 0, 2);
            MPI.COMM_WORLD.Send(eta, 0, eta.length, MPI.DOUBLE, 0, 3);
        }

        x[m] = (eta[m] + epsilon[m] * beta[m - 1]) / (1 - epsilon[m] * alpha[m - 1]);
        x[m - 1] = (beta[m - 1] + alpha[m - 1] * eta[m]) / (1 - epsilon[m] * alpha[m - 1]);

        if (processRank == 0) {
            for (int i = m - 2; i >= 0; i--) {
                x[i] = alpha[i] * x[i + 1] + beta[i];
            }

            MPI.COMM_WORLD.Send(x, 0, x.length, MPI.DOUBLE, 1, 4);
        } else {
            for (int i = m; i <= n - 1 - 1; i++) {
                x[i + 1] = epsilon[i + 1] * x[i] + eta[i + 1];
            }

            double[] xTemp = new double[n];
            MPI.COMM_WORLD.Recv(xTemp, 0, xTemp.length, MPI.DOUBLE, 0, 4);

            for (int i = 0; i <= m - 2; i++) {
                x[i] = xTemp[i];
            }

            MPI.COMM_WORLD.Send(x, 0, x.length, MPI.DOUBLE, 0, 5);
        }
    }

    private static double[] calculateAlpha(double[][] A, int m) {
        double[] vector = new double[m];

        vector[0] = -A[0][1] / A[0][0];

        for (int i = 1; i <= m - 1; i++) {
            vector[i] = -A[i][i + 1] / (A[i][i] + A[i][i - 1] * vector[i - 1]);
        }

        return vector;
    }

    private static double[] calculateBeta(double[][] A, double[] b, int m, double[] alpha) {
        double[] vector = new double[m];

        vector[0] = b[0] / A[0][0];

        for (int i = 1; i <= m - 1; i++) {
            vector[i] = (b[i] - A[i][i - 1] * vector[i - 1]) / (A[i][i] + A[i][i - 1] * alpha[i - 1]);
        }

        return vector;
    }

    private static double[] calculateEpsilon(double[][] A, int n, int m) {
        double[] vector = new double[n];

        vector[n - 1] = -A[n - 1][n - 1 - 1] / A[n - 1][n - 1];

        for (int i = n - 1 - 1; i >= m; i--) {
            vector[i] = -A[i][i - 1] / (A[i][i] + A[i][i + 1] * vector[i + 1]);
        }

        return vector;
    }

    private static double[] calculateEta(double[][] A, double[] b, double[] epsilon, int n, int m) {
        double[] vector = new double[n];

        vector[n - 1] = b[n - 1] / A[n - 1][n - 1];

        for (int i = n - 1 - 1; i >= m; i--) {
            vector[i] = (b[i] - A[i][i + 1] * vector[i + 1]) / (A[i][i] + A[i][i + 1] * epsilon[i + 1]);
        }

        return vector;
    }
}

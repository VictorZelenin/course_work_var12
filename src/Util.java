import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

public class Util {
    public static double[][] decode(double[] array) {
        int size = (int) Math.sqrt(array.length);
        double[][] matrix = new double[size][size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                matrix[i][j] = array[i * size + j];
            }
        }

        return matrix;
    }

    public static double[] encode(double[][] matrix) {
        int size = matrix.length;
        double[] array = new double[size * size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                array[i * size + j] = matrix[i][j];
            }
        }

        return array;
    }

    public static void writeResultInFile(String fileName, double[][] result, double[] vector) {
        File resultFile = new File(fileName);

        if (!resultFile.exists()) {
            try {
                resultFile.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        PrintWriter writer = null;

        try {
            writer = new PrintWriter(resultFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        for (int i = 0; i < vector.length; i++) {
            writer.print(vector[i]);
            writer.write(" ");
            writer.flush();
        }

        for (double[] array : result) {
            writer.write("\r");
            writer.flush();
            for (double d : array) {
                writer.print(d);
                writer.write(" ");
                writer.flush();
            }
        }

        assert writer != null;
        writer.close();
    }

    public static void printVectorInFile(String fileName, double[] vector) {
        File resultFile = new File(fileName);

        if (!resultFile.exists()) {
            try {
                resultFile.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        PrintWriter writer = null;

        try {
            writer = new PrintWriter(resultFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        for (double d : vector) {
            writer.print(d);
            writer.write("\r");
            writer.flush();
        }
    }
}

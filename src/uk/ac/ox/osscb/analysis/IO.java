package uk.ac.ox.osscb.analysis;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

/**
 *
 * @author Michael Golden
 */
public class IO {

    public static void loadFastaSequences(File file, List<String> sequences, List<String> sequenceNames) {
        loadFastaSequences(file, sequences, sequenceNames, Integer.MAX_VALUE);
    }

    public static void loadFastaSequences(File file, List<String> sequences, List<String> sequenceNames, int max) {
        try {
            BufferedReader buffer = new BufferedReader(new FileReader(file));
            String textline = null;
            String sequence = "";
            int n = 0;
            boolean maxReached = false;
            while ((textline = buffer.readLine()) != null) {
                if (maxReached && textline.startsWith(">")) {
                    break;
                }

                if (textline.startsWith(">")) {
                    n++;
                    if (n >= max) {
                        maxReached = true;
                    }

                    sequenceNames.add(textline.substring(1));
                    if (!sequence.equals("")) {
                        sequences.add(sequence.toUpperCase());
                        sequence = "";
                    }
                } else {
                    sequence += textline.trim();
                }

            }
            buffer.close();
            if (!sequence.equals("")) {
                sequences.add(sequence.toUpperCase());
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public static int countFastaSequences(File file) throws IOException {
        int n = 0;
        BufferedReader buffer = new BufferedReader(new FileReader(file));
        String textline = null;
        String sequence = "";

        boolean maxReached = false;
        while ((textline = buffer.readLine()) != null) {
            if (maxReached && textline.startsWith(">")) {
                break;
            }

            if (textline.startsWith(">")) {
                n++;
                /**
                 * if (n >= max) { maxReached = true; }
                 *
                 * sequenceNames.add(textline.substring(1)); if
                 * (!sequence.equals("")) {
                 * sequences.add(sequence.toUpperCase()); sequence = ""; }
                 */
            } else {
                sequence += textline.trim();
            }

        }
        buffer.close();
        /**
         * if (!sequence.equals("")) { sequences.add(sequence); }
         */
        return n;
    }

    public static boolean saveToFASTAfile(List<String> sequences,List<String> sequenceNames, File file) {
        try {
            BufferedWriter buffer = new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < sequences.size() && i < sequenceNames.size(); i++) {
                buffer.write(">" + sequenceNames.get(i) + "\n");
                buffer.write(sequences.get(i) + "\n");
            }
            buffer.close();
            return true;
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return false;
    }

    public static void copyFile(File sourceFile, File destFile) {
        try {
            InputStream in = new FileInputStream(sourceFile);
            OutputStream out = new FileOutputStream(destFile);

            byte[] buf = new byte[1024];
            int len;
            while ((len = in.read(buf)) > 0) {
                out.write(buf, 0, len);
            }
            in.close();
            out.close();
        } catch (FileNotFoundException ex) {
            System.exit(0);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static File getFileInCwd(File cwd, File filename) {
        return new File(cwd.getPath() + "/" + filename.getName());
    }

    public static String[] getHeadersFromCSV(File csvFile) {
        String[] headers = {};
        try {
            BufferedReader buffer = new BufferedReader(new FileReader(csvFile));
            String textline = buffer.readLine();
            if (textline != null) {
                headers = textline.split(",");
            }
            buffer.close();
        } catch (IOException ex) {
            return headers;
        }
        return headers;
    }

    public static ArrayList<String> getColumnFromCSV(File csvFile, int column, boolean hasHeader) {
        ArrayList<String> values = new ArrayList<String>();
        try {
            BufferedReader buffer = new BufferedReader(new FileReader(csvFile));
            String textline = null;
            String[] split = {};
            if (hasHeader) {
                buffer.readLine();
            }
            while ((textline = buffer.readLine()) != null) {
                split = textline.split(",");
                if (column < split.length) {
                    values.add(split[column]);
                } else {
                    values.add(null);
                }
            }
            buffer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return values;
    }

    public static double[] getMinAndMaxValueFromCSV(File csvFile, int column, boolean hasHeader) {
        double[] minAndMax = new double[2];
        minAndMax[0] = Double.MAX_VALUE;
        minAndMax[1] = Double.MIN_VALUE;
        ArrayList<String> data = getColumnFromCSV(csvFile, column, hasHeader);
        for (int i = 0; i < data.size(); i++) {
            try {
                if (data.get(i) != null) {
                    double val = Double.parseDouble(data.get(i));
                    minAndMax[0] = Math.min(minAndMax[0], val);
                    minAndMax[1] = Math.max(minAndMax[1], val);
                }
            } catch (NumberFormatException ex) {
            }
        }
        return minAndMax;
    }

    public static void zipFolder(File inFolder, File outFile, int compressionLevel) throws Exception {
        ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outFile)));
        out.setLevel(Math.min(Math.max(0, compressionLevel), 9));
        int len = inFolder.getAbsolutePath().lastIndexOf(File.separator);
        addFolderToZipStream(inFolder, out, inFolder.getAbsolutePath().substring(0, len + 1));
        out.close();
    }

    private static void addFolderToZipStream(File folder, ZipOutputStream zipStream, String zipFilePath) throws IOException {
        File[] files = folder.listFiles();
        for (int i = 0; i < files.length; i++) {
            if (files[i].isDirectory()) {
                addFolderToZipStream(files[i], zipStream, zipFilePath);
            } else {
                String name = files[i].getAbsolutePath().substring(zipFilePath.length());
                ZipEntry zipEntry = new ZipEntry(name);

                zipStream.putNextEntry(zipEntry);
                BufferedInputStream fileInStream = new BufferedInputStream(new FileInputStream(files[i]));
                BufferedOutputStream zipOutStream = new BufferedOutputStream(zipStream);
                byte[] buffer = new byte[1024];
                while (fileInStream.available() > 0) {
                    fileInStream.read(buffer);
                    zipOutStream.write(buffer);
                }
                fileInStream.close();
                zipStream.closeEntry();
            }
        }
    }

    // http://stackoverflow.com/questions/981578/how-to-unzip-files-recursively-in-java
    public static void extractZip(File zipFile, File extractToPath) throws Exception {
        int BUFFER = 1024;

        ZipFile zip = new ZipFile(zipFile);
        //String newPath = zipFile.substring(0, zipFile.length() - 4);

        extractToPath.mkdirs();
        Enumeration zipFileEntries = zip.entries();

        // Process each entry
        while (zipFileEntries.hasMoreElements()) {
            // grab a zip file entry
            ZipEntry entry = (ZipEntry) zipFileEntries.nextElement();
            String currentEntry = entry.getName();
            File destFile = new File(extractToPath.getAbsolutePath(), currentEntry);
            File destinationParent = destFile.getParentFile();

            // create the parent directory structure if needed
            destinationParent.mkdirs();

            if (!entry.isDirectory()) {
                BufferedInputStream is = new BufferedInputStream(zip.getInputStream(entry));
                int currentByte;
                // establish buffer for writing file
                byte data[] = new byte[BUFFER];

                // write the current file to disk
                FileOutputStream fos = new FileOutputStream(destFile);
                BufferedOutputStream dest = new BufferedOutputStream(fos,
                        BUFFER);

                // read and write until last byte is encountered
                while ((currentByte = is.read(data, 0, BUFFER)) != -1) {
                    dest.write(data, 0, currentByte);
                }
                dest.flush();
                dest.close();
                is.close();
            }

            if (currentEntry.endsWith(".zip")) {
                // found a zip file, try to open
                extractZip(destFile, extractToPath);
            }
        }
    }

    public static boolean contentEquals(File file1, File file2) {
        try {
            BufferedInputStream input1 = new BufferedInputStream(new FileInputStream(file1));
            BufferedInputStream input2 = new BufferedInputStream(new FileInputStream(file2));

            int ch = input1.read();
            while (-1 != ch) {
                int ch2 = input2.read();
                if (ch != ch2) {
                    input1.close();
                    input2.close();
                    return false;
                }
                ch = input1.read();
            }

            int ch2 = input2.read();
            input1.close();
            input2.close();
            return (ch2 == -1);
        } catch (IOException ex) {
            return false;
        }
    }

    public static void main(String[] args) {
        try {
            IO.zipFolder(new File("testcollection2"), new File("testcollection3.zip"), 9);
            IO.extractZip(new File("testcollection3.zip"), new File("C:/myfolder/"));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static void resaveAsNormalisedFasta(File inFastaFile, File normalisedFastaFile) {
        try {
            BufferedReader buffer = new BufferedReader(new FileReader(inFastaFile));
            BufferedWriter writer = new BufferedWriter(new FileWriter(normalisedFastaFile));
            String textline = null;
            String sequence = "";
            int n = 0;
            boolean maxReached = false;
            while ((textline = buffer.readLine()) != null) {
                if (maxReached && textline.startsWith(">")) {
                    break;
                }

                if (textline.startsWith(">")) {     
                    if (!sequence.equals("")) {                        
                        writer.write(">seq"+n+"\n");
                        writer.write(sequence.toUpperCase()+"\n");
                        sequence = "";
                        n++;
                    }
                } else {
                    sequence += textline.trim();
                }

            }
            buffer.close();
            if (!sequence.equals("")) {
                writer.write(">seq"+n+"\n");
                writer.write(sequence.toUpperCase());
                sequence = "";
            }
            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}

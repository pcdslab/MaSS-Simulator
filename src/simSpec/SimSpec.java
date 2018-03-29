/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//package specsimulate;

//import com.sun.org.apache.xpath.internal.operations.Equals;

package simSpec;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.omg.CORBA.CTX_RESTRICT_SCOPE;

import simSpec.MyPair;

/**
 *
 * @author Gul
 */
public class SimSpec {

    /**
     * @param args the command line arguments
     */
    static ArrayList<Double> ionOffsetsC = new ArrayList<>();
    static ArrayList<Double> ionOffsetsN = new ArrayList<>();
    static ArrayList<Double> ionProbsC = new ArrayList<>();
    static ArrayList<Double> ionProbsN = new ArrayList<>();
    static ArrayList<Integer> ionCharges = new ArrayList<>();
    static ArrayList<Double> ionChProbs = new ArrayList<>();
    static ArrayList<Double> ionIntensitiesC = new ArrayList<>();
    static ArrayList<Double> ionIntensitiesN = new ArrayList<>();
    static Double intenCon = 1000.0;
    static Double noiseIntenMin = 0.25;
    static Double noiseIntenMax = 2.25;
    static int totCh = 2;
    static int noiseType = 1; // 0 = no noise, 1 = uniform, 2 = gaussian
    static Double SNR = 4.0;
    static ArrayList<String> immoniumIonsAA = new ArrayList<>();
    static ArrayList<Double> ImmProbs = new ArrayList<>();
    static ArrayList<Double> ImmIntens = new ArrayList<>();
    static Double minMass = 800.0;;
    static Double maxMass = 4000.0;

	static float cTermOff=0, nTermOff=0;
	
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Random r = new Random(System.currentTimeMillis());
        Map<String, Double> aaMasses = new HashMap<>();
        String paramFile;
        String ptmFile;
		//System.out.println("out:::"+args[3]);
        if(Integer.parseInt(args[0])== 0)
            paramFile = "./params.txt";
        else
            paramFile = args[4];
        if(Integer.parseInt(args[1]) == 0)
            ptmFile = "./modifications.ptm";
        else
            ptmFile = args[5];
        String gndTruth = "./peptides.rst";
      //  int myCheck = Integer.parseInt(args[5]);
        BufferedReader ptmReader = new BufferedReader(new FileReader(ptmFile));
        BufferedReader paramReader = new BufferedReader(new FileReader(paramFile));
        readParams(paramReader); // read in the parameter file
       // System.out.println("SNR:"+SNR);
        populateMassTable(aaMasses, ptmReader); // read in PTM table
        String fileName = args[2];//"yeastProteome_digested_Mass400to6000.txt";
        BufferedReader peptideReader = new BufferedReader(new FileReader(fileName));
        String outFile = args[3];
        FileWriter myPepsFile = new FileWriter(outFile);
        PrintWriter printWriter = new PrintWriter(myPepsFile);
        FileWriter myRstFile = new FileWriter(gndTruth);
        PrintWriter rstWriter = new PrintWriter(myRstFile);
        ArrayList<Ion> bList = new ArrayList<>(); // also contains water and ammonia losses
        ArrayList<Ion> yList = new ArrayList<>(); // also contains water and ammonia losses
        ArrayList<Ion> ImmList = new ArrayList<>();
        ArrayList<Ion> noiseList = new ArrayList<>(); // contains noise
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(3);
        df.setGroupingUsed(false);
        
        Random noiseGen = new Random(System.currentTimeMillis());
		Random intenGen = new Random(System.currentTimeMillis());

        String myLine = null;
        boolean checkHeader = false;
        int count = 0;

        printWriter.println("H\tCreationDate\t");
        printWriter.println("H\tSNR\t"+SNR);
        while ((myLine = peptideReader.readLine()) != null) {
            if (checkHeader == false) // if its false then its header
            {
                checkHeader = true;
            } else {
                //String[] myStrArr = myLine.split("\t");
              //  Double peptideMass = Double.parseDouble(myStrArr[3]);
                double peptideMassN = 18.015+cTermOff+nTermOff;
                
                double bSeries = 1+nTermOff;
                //double tempB = 0;
                double ySeries = 0;
               
                int subScriptB = 0;
                int subScriptY = 0;
                int totBions, totYions, ionsGen;
                String peptideStr = myLine;//myStrArr[1];
                
                MyPair c;
                int k = 0;
                while(k < peptideStr.length()){
                	//System.out.println(c);
                	c = getNextAA(peptideStr, k);
                    peptideMassN += aaMasses.get(c.myStr);
                    k = c.myInt + 1;
                }
                
                
                
                if (peptideMassN > minMass && peptideMassN < maxMass) {
					count++;
				//	System.out.println("peptideMass:"+peptideMassN);
				//	System.out.println("scan:"+count);
                    totYions = 0;
                    totBions = 0;
                    ionsGen = 0;
                    ySeries = peptideMassN - bSeries + 2;
                    System.out.println(ySeries);
                    // add ammomnium ions
                    for(int i = 0; i < immoniumIonsAA.size(); i++) {
                    	if(getSelection(ImmProbs.get(i), r));
                    		ImmList.add(new Ion(aaMasses.get(immoniumIonsAA.get(i))-26.99, ImmIntens.get(i)*intenCon));
                    }
                    
                    
                    for (int i = 0; i < ionOffsetsC.size(); i++) {
                        Ion tempY = new Ion(ySeries + ionOffsetsC.get(i), ionIntensitiesC.get(i)*intenCon);
                        //tempY = ySeries + ionOffsetsC.get(i);
                        if (tempY.m_z > 0 && getSelection(ionProbsC.get(i), r)) {
                            yList.add(tempY);
                            ionsGen++;
                            if(ionOffsetsC.get(i) == 0)
                                totYions++;
                        }
                    }

                    //adding multi charged variants
                    for (int i = 0; i < ionCharges.size(); i++) {
                        Ion tempY = new Ion((ySeries + ionCharges.get(i) - 1) / ionCharges.get(i), (ionIntensitiesC.get(0)/2)*intenCon);
                        if (tempY.m_z > 0 && getSelection(ionChProbs.get(i), r)) {
                            ionsGen++;
                            yList.add(tempY);
                        }
                    }
             
                    //*************************************************************//
                    int l = 0;
                    MyPair myChar;
                    while ( l < peptideStr.length()) {
                    	myChar = getNextAA(peptideStr, l);
                        subScriptB++;
                        subScriptY = peptideStr.length() - subScriptB;

                        bSeries += aaMasses.get(myChar.myStr);
                        ySeries = peptideMassN - bSeries + 2;// + cTermOff;

                        for (int i = 0; i < ionOffsetsN.size(); i++) {
                            Ion tempB = new Ion(bSeries + ionOffsetsN.get(i), ionIntensitiesN.get(i)*intenCon);
                            if (tempB.m_z > 0 && getSelection(ionProbsN.get(i), r)) {
                                ionsGen++;
                                bList.add(tempB);
                                if(ionOffsetsN.get(i) == 0)
                                    totBions++;
                            }
                        }
                        //adding multi charged variants
                        for (int i = 0; i < ionCharges.size(); i++) {
                            Ion tempB = new Ion((bSeries + ionCharges.get(i) - 1) / ionCharges.get(i), (ionIntensitiesN.get(0)/2)*intenCon);
                            if (tempB.m_z > 0 && getSelection(ionChProbs.get(i), r)) {
                                bList.add(tempB);
                                ionsGen++;
                            }
                        }

                        for (int i = 0; i < ionOffsetsC.size(); i++) {
                            Ion tempY = new Ion(ySeries + ionOffsetsC.get(i), ionIntensitiesC.get(i)*intenCon);
                            if (tempY.m_z > 0 && getSelection(ionProbsC.get(i), r)) {
                                yList.add(tempY);
                                ionsGen++;
                                if(ionOffsetsC.get(i) == 0)
                                    totYions++;
                                        
                            }
                        }

                        //adding multi charged variants
                        for (int i = 0; i < ionCharges.size(); i++) {
                            Ion tempY = new Ion((ySeries + ionCharges.get(i) - 1) / ionCharges.get(i), (ionIntensitiesC.get(0)/2)*intenCon );
                            if (tempY.m_z > 0 && getSelection(ionChProbs.get(i), r)) {
                                yList.add(tempY);
                                ionsGen++;
                            }
                        }
                        l = myChar.myInt + 1;
                    }

                    //******************************//
                    
                   // System.out.println("IonsGen:"+ionsGen+" B-ions:"+totBions+" Y-ions:"+totYions);
                    int totPeaksReq = (int) ((double)(100*(totBions+totYions))/SNR);
                
                    int totNoisePeaks = totPeaksReq - ionsGen;
                   // System.out.println("totNoise:"+totNoisePeaks+" totPeaksReq:"+totPeaksReq+" ionsGen:"+ionsGen+" totB:"+totBions+" totY:"+totYions);
                    bList.removeAll(yList);
                    bList.removeAll(ImmList);
                    bList.addAll(yList);
                    bList.addAll(ImmList);
                    addNoise(noiseType, noiseList, totNoisePeaks, peptideMassN, noiseGen, intenGen);
                    noiseList.removeAll(bList);
                    noiseList.addAll(bList);
                    noiseList.sort(comp);

                    int maxCharge = (ionCharges.size() > 0)?(ionCharges.get(ionCharges.size() - 1)):1;

                    rstWriter.println("scan:"+count+" peptide:"+peptideStr);
                    printWriter.println("S\t"+count+"\t"+count+"\t"+df.format((peptideMassN+maxCharge)/maxCharge));
                    printWriter.println("Z\t"+maxCharge+"\t"+df.format(peptideMassN+1));
                    for (Ion b : noiseList) {
						if(b.m_z != 0)
							printWriter.println(df.format(b.m_z)+" "+df.format(b.intensity));
                    }
                }
            }
            
            if (count > 50000) {
                break;
            }
			bList.clear();
            yList.clear();
            noiseList.clear();
            ImmList.clear();
        }
        printWriter.flush();
        printWriter.close();
        rstWriter.flush();
        rstWriter.close();
        System.out.println("done");
    }

    static void populateMassTable(Map<String, Double> aaMasses, BufferedReader ptmFile) throws IOException {
        aaMasses.put("G", 57.02);
        aaMasses.put("A", 71.037);
        aaMasses.put("S", 87.03);
        aaMasses.put("P", 97.052);
        aaMasses.put("V", 99.068);
        aaMasses.put("T", 101.047);
        aaMasses.put("C", 103.009);
        aaMasses.put("I", 113.08);
        aaMasses.put("L", 113.08);
        aaMasses.put("N", 114.04);
        aaMasses.put("D", 115.03);
        aaMasses.put("Q", 128.058);
        aaMasses.put("K", 128.09);
        aaMasses.put("E", 129.042);
        aaMasses.put("M", 131.04);
        aaMasses.put("H", 137.06);
        aaMasses.put("F", 147.07);
        aaMasses.put("R", 156.188);
        aaMasses.put("Y", 163.06);
        aaMasses.put("W", 186.07);
        String[] modMass;
        String myLine;
        while((myLine = ptmFile.readLine()) != null){
        	String[] modTok = myLine.split("=");
        	if(modTok[0].equals("cTerm")) {
        		cTermOff = Float.parseFloat(modTok[1]);
        	}else if(modTok[0].equals("nTerm")) {
        		nTermOff = Float.parseFloat(modTok[1]);
        	}
        	else {
            modMass = modTok[1].split("\\+");
            aaMasses.put(modTok[0], aaMasses.get(modMass[0])+Double.parseDouble(modMass[1]));
        	}
        }
    }

    static Comparator<Ion> comp = (Ion a, Ion b) -> {
        return a.compareTo(b);
    };
    
   

    static boolean getSelection(Double prob, Random r) {
        int selection = 0;
        Double threshold = 10 - prob / 10;
        Double randomNum = 0 + r.nextDouble() * 10;
        selection = (randomNum < threshold) ? 0 : 1;

        return selection == 1;
    }

    static void readParams(BufferedReader myFile) throws IOException {
        String myLine;
        String[] lineToks;
        while ((myLine = myFile.readLine()) != null) {
            lineToks = myLine.split(" ");
            if ("C-Ions".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionOffsetsC.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("N-Ions".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionOffsetsN.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("C-Probabilities".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionProbsC.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("N-Probabilities".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionProbsN.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("Max-Charge".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionCharges.add(Integer.parseInt(lineToks[i]));
                }
            } else if ("Prob-Charge".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionChProbs.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("C-Intensities".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionIntensitiesC.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("N-Intensities".equals(lineToks[0])) {
                for (int i = 1; i < lineToks.length; i++) {
                    ionIntensitiesN.add(Double.parseDouble(lineToks[i]));
                }
            } else if ("SNR".equals(lineToks[0])) {
                SNR = Double.parseDouble(lineToks[1]);
            }else if ("Noise".equals(lineToks[0])) {
                noiseType = Integer.parseInt(lineToks[1]);
            }else if ("Noise-Intensities".equals(lineToks[0])) {
                noiseIntenMin = Double.parseDouble(lineToks[1]);
                noiseIntenMax = Double.parseDouble(lineToks[2]);
            }else if("Immonium-Ions".equals(lineToks[0])){
            	for(int i = 1; i < lineToks.length; i++)
            		immoniumIonsAA.add(lineToks[i]);
            }else if("Immonium-Probabilities".equals(lineToks[0])) {
            	for(int i = 1; i < lineToks.length; i++)
            		ImmProbs.add(Double.parseDouble(lineToks[i]));
            }else if("Immonium-Intensities".equals(lineToks[0])) {
            	for(int i = 1; i < lineToks.length; i++)
            		ImmIntens.add(Double.parseDouble(lineToks[i]));
            }else if("Mass-Range".equals(lineToks[0])) {
            		minMass = Double.parseDouble(lineToks[1]);
            		maxMass = Double.parseDouble(lineToks[2]);
            }

        }

        if ((ionOffsetsC.size() != ionProbsC.size()) || (ionOffsetsN.size() != ionProbsN.size()) || (ionChProbs.size() != ionCharges.size())) {
            System.out.println("*****ERROR in params.txt file, please check formatting guide*******");
            System.exit(1);
        }
    }
    
    
    static void addNoise(int noiseType, ArrayList<Ion> noiseList, int totNoisePeaks, Double peptideMass, Random noiseGen, Random intenGen) {
        int noiseAdded = 0;
        Double noiseMean = peptideMass / 2;
        Double noiseStdDev = noiseMean / 2;
		
        if (noiseType == 2) {

            while (noiseAdded < totNoisePeaks) {
                double currentIon = noiseGen.nextGaussian() * noiseStdDev + noiseMean;
				Double noiseInten = noiseIntenMin + intenGen.nextDouble()*(noiseIntenMax-noiseIntenMin);
                if (!noiseList.contains(currentIon) && currentIon > 1 && currentIon < peptideMass + 50) {
                    noiseList.add(new Ion(currentIon, noiseInten*intenCon));
                    noiseAdded++;
                }
            }

        }
        else if(noiseType == 1){
                while (noiseAdded < totNoisePeaks) {
                double currentIon = noiseGen.nextDouble()*peptideMass;
				Double noiseInten = noiseIntenMin + intenGen.nextDouble()*(noiseIntenMax-noiseIntenMin);
                if (!noiseList.contains(currentIon) && currentIon > 1 && currentIon < peptideMass + 50) {
                   noiseList.add(new Ion(currentIon, noiseInten*intenCon));
                    noiseAdded++;
                }
            }
        }
	else if(noiseType == 0){
		// do nothing	
	}
    }
    
	static MyPair getNextAA(String myString, int i) {
		MyPair myPair = new MyPair() ;
		String[] rsult = myString.split("");
		String temp = null;
		String newStr = null;
			temp = rsult[i];
			if(i+1 < rsult.length) {
				if(rsult[i+1].equals("[")) {
				newStr = temp;
				int j = i+1;
						while(!rsult[j].equals("]")) {
							newStr += rsult[j];
							j++;
						}
						newStr += rsult[j];
						temp = newStr;
						myPair.myStr = temp;
						myPair.myInt = j;
						myPair.mod = true;
						//i = j;
			}
		}
			if(myPair.mod == false) {
				myPair.myStr = temp;
				myPair.myInt = i;
			}
		return myPair;
	}
    
    
}

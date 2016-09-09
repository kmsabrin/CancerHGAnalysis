import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;


public class CancerData {
	HashMap<String, ArrayList<Double>> expData429;
	HashMap<String, ArrayList<Double>> expDataNC;
//	HashMap<String, ArrayList<Double>> normalizedExpValues429;
//	HashMap<String, ArrayList<Double>> normalizedExpValuesNC;
//	HashMap<String, ArrayList<Double>> stageValues429;
//	HashMap<String, ArrayList<Double>> stageValuesNC;
	HashMap<String, ArrayList<Double>> transitionDeltas429;
	HashMap<String, ArrayList<Double>> transitionDeltasNC;
	HashMap<String, HashSet<String>> geneProbesetMap;
	
	HashMap<Integer, HashSet<String>> transitoryGenes429;
	HashMap<Integer, HashSet<String>> transitoryGenesNC;
	
	int[] hourStartIndices = {0, 3, 6, 9, 12, 15};
	int nProbeSets;
	int nReplicas = 3;
	int nStages = 6;
	int nTransitions = 5;
	
	ArrayList<Double> ambientNoiseDistribution;
	ArrayList<Double> transitionDeltaDistribution;
	
	public CancerData() {
		expData429 = new HashMap();
		expDataNC = new HashMap();
		geneProbesetMap = new HashMap();
		
		ambientNoiseDistribution = new ArrayList();
		
		transitionDeltaDistribution = new ArrayList();
		transitionDeltas429 = new HashMap();
		transitionDeltasNC = new HashMap();
		
		transitoryGenes429 = new HashMap();
		transitoryGenesNC = new HashMap();
	}
	
	public void insertMappedList(HashMap<String, ArrayList<Double>> mappedList, String key, Double value) {
		if (mappedList.containsKey(key)) {
			mappedList.get(key).add(value);
			return;
		}
		ArrayList<Double> newList = new ArrayList();
		newList.add(value);
		mappedList.put(key, newList);
	}
	
	public void loadCancerData() throws Exception {
		Scanner scanner = new Scanner(new File("CancerData//entire_dataset.txt"));
		
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			
			String probeSetId = tokens[0];
			String geneName = tokens[1];
			for (int i = 2; i <= 4; ++i) {
				insertMappedList(expData429, probeSetId, Double.parseDouble(tokens[i]));
				insertMappedList(expDataNC, probeSetId, Double.parseDouble(tokens[i]));
			}
			
			for (int i = 5; i <= 19; ++i) {
				insertMappedList(expData429, probeSetId, Double.parseDouble(tokens[i]));
			}
			
			for (int i = 20; i <= 34; ++i) {
				insertMappedList(expDataNC, probeSetId, Double.parseDouble(tokens[i]));
			}
			//special correction for nc-3h-0
			expDataNC.get(probeSetId).set(3, (expDataNC.get(probeSetId).get(4) + expDataNC.get(probeSetId).get(5)) * 0.5);
			
			if (geneProbesetMap.containsKey(geneName)) {
				geneProbesetMap.get(geneName).add(probeSetId);
			}
			else {
				HashSet<String> hset = new HashSet();
				hset.add(probeSetId);
				geneProbesetMap.put(geneName, hset);
			}
			
			++nProbeSets;
		}
		
		scanner.close();
	}
	
	public void checkOffData(HashMap<String, ArrayList<Double>> expValues) {
		for (int hourIndex: hourStartIndices) {
			double leftGapValues[] = new double[nProbeSets];
			double rightGapValues[] = new double[nProbeSets];
			int idx1 = 0;
			for (String probeSetId: expValues.keySet()) {
				double replicatedExpValues[] = new double[nReplicas];
				int idx2 = 0;
				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
					replicatedExpValues[idx2++] = d;
				}
				// to sort or not to sort
				//Arrays.sort(replicatedExpValues);
				
				leftGapValues[idx1] = (replicatedExpValues[1] - replicatedExpValues[0]);// / replicatedExpValues[1];
				rightGapValues[idx1] = (replicatedExpValues[2] - replicatedExpValues[1]);// / replicatedExpValues[1];
				//System.out.println(probeSetId);
				//System.out.println(replicatedExpValues[0] + "\t" + replicatedExpValues[1] + "\t" + replicatedExpValues[2]);
				//System.out.println(leftGapValues[idx1] + "\t" + rightGapValues[idx1]);
				idx1++;
			}
			
			System.out.println(StatUtils.percentile(leftGapValues, 25) 
								+ "\t" + StatUtils.percentile(leftGapValues, 50) 
								+ "\t" + StatUtils.percentile(leftGapValues, 75) 
								+ "\n" + StatUtils.percentile(rightGapValues, 25) 
								+ "\t" + StatUtils.percentile(rightGapValues, 50) 
								+ "\t" + StatUtils.percentile(rightGapValues, 75));
		}
	}
	
	public double getSpearmanRhoTransitionDelta(int sz, String geneName, HashMap<String, ArrayList<Double>> expValues) {
		double transitionDelta[][] = new double[sz][nTransitions];
		int t = 0;
		for (String probeSetId: geneProbesetMap.get(geneName)) {
			ArrayList<Double> expList = expValues.get(probeSetId);
			int k = 0;
			for (int i = 0; i < nTransitions; i += nReplicas) {
				int j = i + nReplicas;
				double expV0 = (expList.get(i) + expList.get(i + 1) + expList.get(i + 2)) / 3.0;
				double expV1 = (expList.get(j) + expList.get(j + 1) + expList.get(j + 2)) / 3.0;
				transitionDelta[t][k++] = expV1 - expV0;		
			}
			t++;
		}
		
		double avgSpearmanRho = 0;
		for (int i = 0; i < sz - 1; ++i) {
			for (int j = i + 1; j < sz; ++j) {
				double spearmanRho = new SpearmansCorrelation().correlation(transitionDelta[i], transitionDelta[j]);
				avgSpearmanRho += spearmanRho;
			}
		}
		avgSpearmanRho /= (sz * (sz - 1) / 2.0);
		return avgSpearmanRho;
	}
	
	public void checkProbesetConsistency() {
		for (String geneName: geneProbesetMap.keySet()) {
			int sz = geneProbesetMap.get(geneName).size(); 
			if (sz < 2) {
				continue;
			}
			double spearmanR429 = getSpearmanRhoTransitionDelta(sz, geneName, expData429);
			double spearmanRNC = getSpearmanRhoTransitionDelta(sz, geneName, expDataNC);
			System.out.println(spearmanR429 + "\t" + spearmanRNC);
		}
	}
	
	public void checkDataConsistency() {
		checkOffData(expData429);
		System.out.println("\n\n");
		checkOffData(expDataNC);
		
//		checkProbesetConsistency();
	}
	
	public void normalizeDataHelper(HashMap<String, ArrayList<Double>> expValues) {
		double normalizer = 0;
		for (int hourIndex: hourStartIndices) {
			for (String probeSetId: expValues.keySet()) {
				double replicatedExpValues[] = new double[nReplicas];
				int idx = 0;
				double avg = 0;
				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
					replicatedExpValues[idx++] = d;
					avg += d;
				}
				avg /= nReplicas;
				normalizer += avg;
			}
			normalizer /= expValues.keySet().size();
			for (String probeSetId: expValues.keySet()) {
				double normalizedValue = expValues.get(probeSetId).get(hourIndex) / normalizer;
				expValues.get(probeSetId).set(hourIndex, normalizedValue);
			}
		}
	}
	
	public void normalizeData() {
		normalizeDataHelper(expData429);
		normalizeDataHelper(expDataNC);
		nReplicas = 1;
	}
	
	public double ambientNoisesHelper(HashMap<String, ArrayList<Double>> expValues) {
		double sum = 0;
		for (int hourIndex: hourStartIndices) {
			for (String probeSetId: expValues.keySet()) {
				double replicatedExpValues[] = new double[nReplicas];
				int idx = 0;
				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
					replicatedExpValues[idx++] = d;
				}
//				if (Math.abs(replicatedExpValues[1] - replicatedExpValues[2]) > 6) {
//					System.out.println(probeSetId + "\t" + hourIndex);
//				}
				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[1] - replicatedExpValues[0]));
				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[2] - replicatedExpValues[0]));
				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[2] - replicatedExpValues[1]));
				sum += Math.abs(replicatedExpValues[1] - replicatedExpValues[0]);
				sum += Math.abs(replicatedExpValues[2] - replicatedExpValues[0]);
				sum += Math.abs(replicatedExpValues[2] - replicatedExpValues[1]);
			}
		}
		return sum;
	}
	
	public void getAmbientNoises() throws Exception {
		double normalizerSum = 0;
		normalizerSum += ambientNoisesHelper(expData429);
		normalizerSum += ambientNoisesHelper(expDataNC);
		
		PrintWriter pw = new PrintWriter(new File("ambientNoiseDistribution.txt"));
		double noiseDistributionArray[] = new double[ambientNoiseDistribution.size()];
		for (int idx = 0; idx < ambientNoiseDistribution.size(); ++idx) {
//			ambientNoiseValues.set(idx, ambientNoiseValues.get(idx) / normalizerSum);
			noiseDistributionArray[idx] = ambientNoiseDistribution.get(idx);
//			System.out.println(noiseDistributionArray[idx]);
			pw.println(noiseDistributionArray[idx]);
		}
		pw.close();
		System.out.println(noiseDistributionArray.length);
		System.out.println(StatUtils.mean(noiseDistributionArray) + "\t" + Math.sqrt(StatUtils.variance(noiseDistributionArray)));
//		System.out.println(StatUtils.percentile(noiseDistributionArray, 10));
//		System.out.println(StatUtils.percentile(noiseDistributionArray, 50));
		System.out.println(StatUtils.percentile(noiseDistributionArray, 95));
	}
	
//	public void getStageValues() {
//		
//	}
//	
//	public void getStageValuesHelper(HashMap<String, ArrayList<Double>> expValues) {
//		for (int hourIndex: hourStartIndices) {
//			for (String probeSetId: expValues.keySet()) {
//				double replicatedExpValues[] = new double[nReplicas];
//				int idx = 0;
//				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
//					replicatedExpValues[idx++] = d;
//				}
//			}
//		}
//	}
	
	public void getTransitionDeltas() throws Exception {
		getTransitionDeltasHelper(expData429, transitionDeltas429);
		getTransitionDeltasHelper(expDataNC, transitionDeltasNC);
		
		PrintWriter pw = new PrintWriter(new File("deltaDistributionNorm.txt"));
		double deltaDistributionArray[] = new double[transitionDeltaDistribution.size()];
		for (int idx = 0; idx < transitionDeltaDistribution.size(); ++idx) {
			deltaDistributionArray[idx] = transitionDeltaDistribution.get(idx);
//			System.out.println(deltaDistributionArray[idx]);
			pw.println(deltaDistributionArray[idx]);
		}
		pw.close();
		System.out.println(deltaDistributionArray.length);
		System.out.println(StatUtils.mean(deltaDistributionArray) + "\t" + Math.sqrt(StatUtils.variance(deltaDistributionArray)));
		System.out.println(StatUtils.percentile(deltaDistributionArray, 95));
		
//		transitionDeltaDistribution.clear();
//		getTransitionDeltasHelper(expDataNC, transitionDeltasNC);
//		deltaDistributionArray = new double[transitionDeltaDistribution.size()];
//		for (int idx = 0; idx < transitionDeltaDistribution.size(); ++idx) {
//			deltaDistributionArray[idx] = transitionDeltaDistribution.get(idx);
//		}
//		System.out.println(deltaDistributionArray.length);
//		System.out.println(StatUtils.mean(deltaDistributionArray) + "\t" + Math.sqrt(StatUtils.variance(deltaDistributionArray)));
	}
	
	public void getTransitionDeltasHelper(HashMap<String, ArrayList<Double>> expValues, HashMap<String, ArrayList<Double>> deltaValues) {
		for (String probeSetId: expValues.keySet()) {
			deltaValues.put(probeSetId, new ArrayList<Double>());
			for (int idx = 0; idx < hourStartIndices.length - 1; ++idx) {
				int h0 = hourStartIndices[idx];
				double replicatedExpValuesH0[] = new double[nReplicas];
				int idx2 = 0;
				double avgH0 = 0;
				for(double d: expValues.get(probeSetId).subList(h0, h0 + nReplicas)) {
					replicatedExpValuesH0[idx2++] = d;
					avgH0 += d;
//					System.out.println(d);
				}
				avgH0 /= nReplicas;
//				System.out.println("avg: " + avgH0);
				
				int h1 = hourStartIndices[idx + 1];
				double replicatedExpValuesH1[] = new double[nReplicas];
				idx2 = 0;
				double avgH1 = 0;
				for(double d: expValues.get(probeSetId).subList(h1, h1 + nReplicas)) {
					replicatedExpValuesH1[idx2++] = d;
					avgH1 += d;
//					System.out.println(d);
				}
				avgH1 /= nReplicas;
//				System.out.println("--- --- ---");
				
				double deltaArray[] = new double[nReplicas * nReplicas];
				int k = 0;
				for (int i = 0; i < nReplicas; ++i) {
					for (int j = 0; j < nReplicas; ++j) {
						double delta = Math.abs(replicatedExpValuesH0[i] - replicatedExpValuesH1[j]);
						deltaArray[k++] = delta;
						transitionDeltaDistribution.add(delta);
//						transitionDeltaDistribution.add(Math.abs(avgH0 - avgH1));
					}
				}
				deltaValues.get(probeSetId).add(StatUtils.percentile(deltaArray, 50));
//				deltaValues.get(probeSetId).add(Math.abs(avgH0 - avgH1));
//				System.out.println(probeSetId + "\t" + h0 + "\t" + Math.abs(avgH0 - avgH1));
//				System.out.println(" --- --- --- ");
			}
		}
	}
	
	public void getTransitoryGenes() {
		for (int i = 0; i < nTransitions; ++i) {
			transitoryGenes429.put(i, new HashSet<String>());
			transitoryGenesNC.put(i, new HashSet<String>());
		}
		
		double threshold = 2.85;
		for (String probeSetId: transitionDeltas429.keySet()) {
			ArrayList<Double> deltas429 = transitionDeltas429.get(probeSetId);
			ArrayList<Double> deltasNC = transitionDeltasNC.get(probeSetId);
			for (int i = 0; i < nTransitions; ++i) {
				if (deltas429.get(i) > threshold) {
					transitoryGenes429.get(i).add(probeSetId);
				}
				if (deltasNC.get(i) > threshold) {
					transitoryGenesNC.get(i).add(probeSetId);
				}
			}
		}
		
		for (int i = 0; i < nTransitions; ++i) {
			System.out.println(transitoryGenes429.get(i).size() + "\t" + transitoryGenesNC.get(i).size());
			
			if (transitoryGenes429.get(i).size() < 30) {
				for (String s: transitoryGenes429.get(i)) {
					System.out.println(s);
				}
				System.out.println("## ## ##");
			}
			if (transitoryGenesNC.get(i).size() < 30) {
				for (String s: transitoryGenesNC.get(i)) {
					System.out.println(s);
				}
				System.out.println("** ** **");
			}
		}
	}
	
	public static void main(String[] args) throws Exception {
		CancerData cancerData = new CancerData();
		cancerData.loadCancerData();
//		cancerData.checkDataConsistency();
		cancerData.getAmbientNoises();
		cancerData.normalizeData();
		cancerData.getTransitionDeltas();
		cancerData.getTransitoryGenes();
	}
}

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.TreeMap;

import org.apache.commons.math3.stat.StatUtils;


public class CancerData {
	int nProbeSets;
	int nReplicas = 9;
	int nStages = 6;
	int nTransitions = 5;
	double transitionThreshold = 2.5;
	int nProbesets = 41580;

	HashMap<String, ArrayList<Double>> expData429;
	HashMap<String, ArrayList<Double>> expDataNC;
	HashMap<String, ArrayList<Double>> ncN429;

	HashMap<String, ArrayList<Double>> transitionDeltas429;
	HashMap<String, ArrayList<Double>> transitionDeltasNC;
	HashMap<String, ArrayList<Double>> transitionDeltasNCn429;
	
	HashMap<String, HashSet<String>> geneProbesetMap;
	HashMap<String, String> probesetGeneMap;
	
	HashMap<Integer, HashSet<String>> transitoryGenes429;
	HashMap<Integer, HashSet<String>> transitoryGenesNC;
	HashMap<Integer, HashSet<String>> transitoryGenesNCn429;
	
	double stageNoiseStD429[] = new double[nStages];
	double stageNoiseStDNC[] = new double[nStages];
	double stageNoiseStDNCn429[] = new double[nStages];
	
	ArrayList<Double> ambientNoiseDistribution;
	ArrayList<Double> transitionDeltaDistribution;
		
	ArrayList<String> emtGenes;
	HashMap<String, Integer> emtGeneTransition429;
	HashMap<String, Integer> emtGeneTransitionNC;
	HashMap<String, Integer> emtGeneTransitionNCn429;
	
	HashMap<Integer, HashSet<String>> significantStageGenes;
	
		
	private void getNCn429() {
		int oReplica = 3;
		for (String s: expData429.keySet()) {
			ncN429.put(s, new ArrayList<Double>());
			for (int stage = 0; stage < nStages; ++stage) {
				List<Double> list429 = expData429.get(s).subList(stage * oReplica, (stage + 1) * oReplica);
				List<Double> listNC = expDataNC.get(s).subList(stage * oReplica, (stage + 1) * oReplica);
				for (double d429: list429) {
					for (double dNC: listNC) {
						if (stage > 0) {
							ncN429.get(s).add(d429 - dNC);
						}
						else {
							ncN429.get(s).add(d429);
						}
//						if (s.equals("1053_at")) {
//							System.out.println(stage + " : " + (d429 - dNC) + " # " + d429 + " # " + dNC );
//						}
					}
				}
			}
		}
	}
	
	public CancerData() {
		expData429 = new HashMap();
		expDataNC = new HashMap();
		ncN429 = new HashMap();
		
		geneProbesetMap = new HashMap();
		probesetGeneMap = new HashMap();
		
		ambientNoiseDistribution = new ArrayList();
		
		transitionDeltaDistribution = new ArrayList();
		transitionDeltas429 = new HashMap();
		transitionDeltasNC = new HashMap();
		transitionDeltasNCn429 = new HashMap();
		
		transitoryGenes429 = new HashMap();
		transitoryGenesNC = new HashMap();
		transitoryGenesNCn429 = new HashMap();
		
		emtGenes = new ArrayList();
		emtGeneTransition429 = new HashMap();
		emtGeneTransitionNC = new HashMap();
		emtGeneTransitionNCn429 = new HashMap();
		
		significantStageGenes = new HashMap();
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
	
	private void getCDFExpValHelper(HashMap<String, ArrayList<Double>> expValues, TreeMap<Double, Double> histogram) {
		for (String probeSetId: expValues.keySet()) {
			for (double d: expValues.get(probeSetId)) {
				if (histogram.containsKey(d)) {
					double v = histogram.get(d);
					histogram.put(d, v + 1);
				}
				else {
					histogram.put(d, 1.0);
				}
			}
		}
	}
	
	public void getCDFExpVal() throws Exception {
		TreeMap<Double, Double> histogram = new TreeMap();
		getCDFExpValHelper(expData429, histogram);
		getCDFExpValHelper(expDataNC, histogram);
		
		double cumSum = 0;
		PrintWriter pw = new PrintWriter(new File("expCDF.txt"));
		for (double d: histogram.keySet()) {
			cumSum += histogram.get(d);
			pw.println(d + "\t" + (cumSum / (expData429.size() * 18 * 2)));
		}
		pw.close();
	}
	
	private void getCDFDeltasHelper(int stage,
			HashMap<String, ArrayList<Double>> deltaValues,
			String type) throws Exception {
		TreeMap<Double, Double> histogram = new TreeMap();
		for (String probeSetId : deltaValues.keySet()) {
			double d = (deltaValues.get(probeSetId).get(stage));
			if (histogram.containsKey(d)) {
				double v = histogram.get(d);
				histogram.put(d, v + 1);
			} else {
				histogram.put(d, 1.0);
			}
		}
		
		double cumSum = 0;
		PrintWriter pw = new PrintWriter(new File("deltaCDF-" + stage + "-" +  type + ".txt"));
		for (double d: histogram.keySet()) {
			cumSum += histogram.get(d);
			pw.println(d + "\t" + (cumSum / deltaValues.size()));
		}
		pw.close();
	}
	
	
	public void getCDFDeltas() throws Exception {
		for (int i = 0; i < nTransitions; ++i) {
//			getCDFDeltasHelper(i, transitionDeltas429, "m429");
//			getCDFDeltasHelper(i, transitionDeltasNC, "NC");
			getCDFDeltasHelper(i, transitionDeltasNCn429, "NCn429");
			
		}
	}
	
	private void getCDFNoiseDistributionHelper(int stage,
			double dist[],
			String type) throws Exception {
		TreeMap<Double, Double> histogram = new TreeMap();
		for (double d : dist) {
			if (histogram.containsKey(d)) {
				double v = histogram.get(d);
				histogram.put(d, v + 1);
			} else {
				histogram.put(d, 1.0);
			}
		}
		
		double cumSum = 0;
		PrintWriter pw = new PrintWriter(new File("noiseCDF-" + stage + "-" +  type + ".txt"));
		for (double d: histogram.keySet()) {
			cumSum += histogram.get(d);
			pw.println(d + "\t" + (cumSum / dist.length));
		}
		pw.close();
	}
	
	private void loadEMTGenes() throws Exception {
		Scanner scanner = new Scanner(new File("CancerData//EMT_genes_X.txt"));
		while (scanner.hasNextLine()) {
			String gene = scanner.nextLine();
			emtGenes.add(gene);
			if (!geneProbesetMap.containsKey(gene)) {
//				System.out.println(gene);
//				findGenes(gene);
			}
		}
		scanner.close();
	}
	
	public void loadCancerData() throws Exception {
		Scanner scanner = new Scanner(new File("CancerData//entire_dataset.txt"));
		
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			
			String probeSetId = tokens[0];
			String gene = tokens[1];
			probesetGeneMap.put(probeSetId, gene);
			
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
			
			// special correction for nc-3h-0
			expDataNC.get(probeSetId).set(3, (expDataNC.get(probeSetId).get(4) + expDataNC.get(probeSetId).get(5)) * 0.5);
			
			if (geneProbesetMap.containsKey(gene)) {
				geneProbesetMap.get(gene).add(probeSetId);
			}
			else {
				HashSet<String> hset = new HashSet();
				hset.add(probeSetId);
				geneProbesetMap.put(gene, hset);
			}
			
			++nProbeSets;
		}
		
//		System.out.println(geneProbesetMap.size());
		scanner.close();	
		
		loadEMTGenes();
		
		getNCn429();
	}
	
	private void checkOffData(HashMap<String, ArrayList<Double>> expValues) {
		for (int stage = 0; stage <nStages; ++stage) {
			int hourIndex = stage * nReplicas;
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
	
	private double getConflictRatioTransitionDelta(int sz, String geneName, HashMap<String, ArrayList<Double>> deltaValues) {
		int transitionDelta[][] = new int[sz][nTransitions];
		int t = 0;
		for (String probeSetId: geneProbesetMap.get(geneName)) {
			ArrayList<Double> deltaList = deltaValues.get(probeSetId);
			int k = 0;
			for (int i = 0; i < nTransitions; i++) {
				int trend = 0;
				if (deltaList.get(i) > transitionThreshold) {
					trend = -1;
				}
				if (-1 * deltaList.get(i) > transitionThreshold) {
					trend = 1;
				}
				transitionDelta[t][k++] = trend;
			}
			t++;
		}
		
		double conflictKount = 0;
		for (int i = 0; i < sz - 1; ++i) {
			for (int j = i + 1; j < sz; ++j) {
				for (int k = 0; k < nTransitions; ++k) {
					if (transitionDelta[i][k] * transitionDelta[j][k] == -1) {
						++conflictKount;
					}
				}
			}
		}
		
		return conflictKount / (sz * (sz - 1.0) * 0.5 * nTransitions);
	}
	
	public void checkProbesetConsistency() {
		System.out.println(geneProbesetMap.size());
		double konflict429 = 0;
		double konflictNC = 0;
		double multiprobeGenes = 0;
		for (String geneName: geneProbesetMap.keySet()) {
			int sz = geneProbesetMap.get(geneName).size(); 
			if (sz < 2) {
				continue;
			}
			++multiprobeGenes;
			double conflictRatio429 = getConflictRatioTransitionDelta(sz, geneName, transitionDeltas429);
			double conflictRatioNC = getConflictRatioTransitionDelta(sz, geneName, transitionDeltasNC);
			if (conflictRatio429 > 0.2) ++konflict429;
			if (conflictRatioNC > 0.2) ++konflictNC;
		}
		
		System.out.println(multiprobeGenes);
		System.out.println((konflict429 / multiprobeGenes) + "\t" + (konflictNC / multiprobeGenes));
	}
	
	public void checkMicroarrayReplicaConsistency() {
		checkOffData(expData429);
		System.out.println("\n\n");
		checkOffData(expDataNC);
	}
	
	private void normalizeDataHelper(HashMap<String, ArrayList<Double>> expValues) {
		for (int stage = 0; stage <nStages; ++stage) {
			int hourIndex = stage * nReplicas;			
			for (int replicaIndex = 0; replicaIndex < nReplicas; ++replicaIndex) {
				double cumSum = 0;
				for (String probeSetId: expValues.keySet()) {
					double d =  expValues.get(probeSetId).get(hourIndex + replicaIndex);
					cumSum += d;
				}
			
				for (String probeSetId: expValues.keySet()) {
					double normalizedValue = expValues.get(probeSetId).get(hourIndex + replicaIndex) / cumSum;
					expValues.get(probeSetId).set(hourIndex + replicaIndex, normalizedValue);
				}
			}
		}
	}
	
	public void normalizeData() {
		normalizeDataHelper(expData429);
		normalizeDataHelper(expDataNC);
	}
	
	private void ambientNoisesHelper(HashMap<String, ArrayList<Double>> expValues) {
		for (int stage = 0; stage <nStages; ++stage) {
			int hourIndex = stage * nReplicas;
			for (String probeSetId: expValues.keySet()) {
				double replicatedExpValues[] = new double[nReplicas];
				int idx = 0;
				double avg = 0;
				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
					replicatedExpValues[idx++] = d;
					avg += d;
				}
				avg /= nReplicas;
//				if (Math.abs(replicatedExpValues[1] - replicatedExpValues[2]) > 6) {
//					System.out.println(probeSetId + "\t" + hourIndex);
//				}
				
//				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[1] - replicatedExpValues[0]));
//				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[2] - replicatedExpValues[0]));
//				ambientNoiseDistribution.add(Math.abs(replicatedExpValues[2] - replicatedExpValues[1]));
				
				ambientNoiseDistribution.add(avg - replicatedExpValues[0]);
				ambientNoiseDistribution.add(avg - replicatedExpValues[1]);
				ambientNoiseDistribution.add(avg - replicatedExpValues[2]);
			}
		}
	}
	
	private double stageNoiseDistributionHelper(int stage, HashMap<String, ArrayList<Double>> expValues, String id) throws Exception {
		double noiseDistribution[] = new double[expValues.size() * nReplicas];
		int k = 0;
		for (String probeSetId: expValues.keySet()) {
			double replicatedExpValues[] = new double[nReplicas];
			int idx = 0;
			double avg = 0;
			for (double d: expValues.get(probeSetId).subList(stage * nReplicas, (stage + 1) * nReplicas)) {
				replicatedExpValues[idx++] = d;
				avg += d;
			}
			avg /= nReplicas;
			for (int i = 0; i < nReplicas; ++i) {
				noiseDistribution[k++] = avg - replicatedExpValues[i];
			}
		}
		
//		getCDFNoiseDistributionHelper(stage, noiseDistribution, id);
		return Math.sqrt(StatUtils.variance(noiseDistribution));
	}
	
	public void getStageNoiseDistribution() throws Exception {
		for (int i = 0; i < nStages; ++i) {
			if (nReplicas > 3) {
				stageNoiseStDNCn429[i] =  stageNoiseDistributionHelper(i, ncN429, "NCn429");
//				System.out.println("(S)igma stage: " + i + " ncN429: " + stageNoiseStDNCn429[i]);
			}
			else {
				stageNoiseStD429[i] =  stageNoiseDistributionHelper(i, expData429, "m429");
				stageNoiseStDNC[i] =  stageNoiseDistributionHelper(i, expDataNC, "NC");
//				System.out.println("(S)igma stage: " + i + " m429: " + stageNoiseStD429[i]  + " nC: " + stageNoiseStDNC[i]);
			}
		}
	}
	
	public void getAmbientNoises() throws Exception {
		ambientNoisesHelper(expData429);
		ambientNoisesHelper(expDataNC);
		
		PrintWriter pw = new PrintWriter(new File("ambientNoiseDistribution.txt"));
		double noiseDistributionArray[] = new double[ambientNoiseDistribution.size()];
		for (int idx = 0; idx < ambientNoiseDistribution.size(); ++idx) {
			noiseDistributionArray[idx] = ambientNoiseDistribution.get(idx);
//			System.out.println(noiseDistributionArray[idx]);
			pw.println(noiseDistributionArray[idx]);
		}
		pw.close();
		System.out.println(noiseDistributionArray.length);
		System.out.println(StatUtils.mean(noiseDistributionArray) + "\t" + 4 * Math.sqrt(StatUtils.variance(noiseDistributionArray)));
//		System.out.println(StatUtils.percentile(noiseDistributionArray, 10));
//		System.out.println(StatUtils.percentile(noiseDistributionArray, 50));
//		System.out.println(StatUtils.percentile(noiseDistributionArray, 95));
	}
	
	public void getTransitionDeltas() throws Exception {
		if (nReplicas > 3) {
			getTransitionDeltasHelper(ncN429, transitionDeltasNCn429);
		}
		else {
			getTransitionDeltasHelper(expData429, transitionDeltas429);
			getTransitionDeltasHelper(expDataNC, transitionDeltasNC);
		}	

//		double deltaDistributionArray[] = new double[transitionDeltaDistribution.size()];
//		for (int idx = 0; idx < transitionDeltaDistribution.size(); ++idx) {
//			deltaDistributionArray[idx] = transitionDeltaDistribution.get(idx);
//		}
//		System.out.println(deltaDistributionArray.length);
//		System.out.println(StatUtils.mean(deltaDistributionArray) + "\t" + Math.sqrt(StatUtils.variance(deltaDistributionArray)));
//		System.out.println(StatUtils.percentile(deltaDistributionArray, 50));

		
//		PrintWriter pw = new PrintWriter(new File("deltaDistributionNorm.txt"));
//		double deltaDistributionArray[] = new double[transitionDeltaDistribution.size()];
//		for (int idx = 0; idx < transitionDeltaDistribution.size(); ++idx) {
//			deltaDistributionArray[idx] = transitionDeltaDistribution.get(idx);
////			System.out.println(deltaDistributionArray[idx]);
//			pw.println(deltaDistributionArray[idx]);
//		}
//		pw.close();
		
//		System.out.println(deltaDistributionArray.length);
//		System.out.println(StatUtils.mean(deltaDistributionArray) + "\t" + Math.sqrt(StatUtils.variance(deltaDistributionArray)));
//		System.out.println(StatUtils.percentile(deltaDistributionArray, 95));
		
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
			for (int stage = 0; stage < nStages - 1; ++stage) {
				int h0 = stage * nReplicas;
				double replicatedExpValuesH0[] = new double[nReplicas];
				int idx = 0;
				for(double d: expValues.get(probeSetId).subList(h0, h0 + nReplicas)) {
					replicatedExpValuesH0[idx++] = d;
				}
				
				int h1 = (stage + 1) * nReplicas;
				double replicatedExpValuesH1[] = new double[nReplicas];
				idx = 0;
				for(double d: expValues.get(probeSetId).subList(h1, h1 + nReplicas)) {
					replicatedExpValuesH1[idx++] = d;
				}
				
				double deltaArray[] = new double[nReplicas * nReplicas];
				int k = 0;
				for (int i = 0; i < nReplicas; ++i) {
					for (int j = 0; j < nReplicas; ++j) {
						double delta = Math.abs(replicatedExpValuesH1[i] - replicatedExpValuesH0[j]);
						deltaArray[k++] = delta;
//						transitionDeltaDistribution.add(delta);
					}
				}

				if (stage > 0) {
//					deltaArray = Arrays.copyOfRange(deltaArray, 0, k); // for nc-3a // fixed
					deltaValues.get(probeSetId).add(StatUtils.percentile(deltaArray, 50));
//					deltaValues.get(probeSetId).add(StatUtils.mean(deltaArray));
				}
				else { // special case for first transition
					for (int i = 0; i < replicatedExpValuesH1.length; ++i) {
						replicatedExpValuesH1[i] = Math.abs(replicatedExpValuesH1[i]);
					}
					deltaValues.get(probeSetId).add(StatUtils.percentile(replicatedExpValuesH1, 50));
				}
			}
		}
	}
	
	public void getSignificantStageValues(
			HashMap<String, ArrayList<Double>> expValues) {
		for (int stage = 1; stage < nStages; ++stage) {
//			System.out.println("Stage " + stage);
			HashMap<String, Double> stageValueMap = new HashMap();
			for (String probeSetId : expValues.keySet()) {
				int h0 = stage * nReplicas;
				double replicatedExpValuesH0[] = new double[nReplicas];
				int idx = 0;
				for (double d : expValues.get(probeSetId).subList(h0, h0 + nReplicas)) {
					replicatedExpValuesH0[idx++] = Math.abs(d);
				}
				double median = StatUtils.percentile(replicatedExpValuesH0, 50);
				stageValueMap.put(probeSetId, median);
			}
			
			double w = 3;
			double threshold = Math.sqrt(stageNoiseStD429[stage] * stageNoiseStD429[stage] + stageNoiseStDNC[stage] * stageNoiseStDNC[stage]);
			threshold *= w;
			HashSet<String> emtSet = new HashSet();
			HashSet<String> transGeneSet = new HashSet();
			HashSet<String> downEMTGeneSet = new HashSet();
			for (String probeSetId : stageValueMap.keySet()) {
				double v = stageValueMap.get(probeSetId);
				if (Math.abs(v) < threshold) {
					continue;
				}
				transGeneSet.add(probesetGeneMap.get(probeSetId));
				if (emtGenes.contains(probesetGeneMap.get(probeSetId))) {
					emtSet.add(probesetGeneMap.get(probeSetId));
//					System.out.println(probesetGeneMap.get(probeSetId));
					if (v < 0) {
						downEMTGeneSet.add(probesetGeneMap.get(probeSetId));
					}
				}
			}
//			System.out.println(
////					stage 
//					(transGeneSet.size() * 1.0 / geneProbesetMap.size()) 
//					+ "\t" + transGeneSet.size()
//					+ "\t" + (emtSet.size() * 1.0 / emtGenes.size())
//					+ "\t" + emtSet.size()
//					+ "\t" 
//					+ (downEMTGeneSet.size() * 1.0 / emtSet.size()));
			
//			for (String s: transGeneSet) {
//				System.out.println(s);
//			}
//			
//			System.out.println("### ### ###");
			
			significantStageGenes.put(stage, new HashSet(transGeneSet));
		}
	}
	
	private int getUpDownRegulation(int stage, 
			                        HashMap<Integer, HashSet<String>> transitoryGenes, 
			                        HashMap<String, ArrayList<Double>> transitionDeltas) {
		// get up/down regulation info
		int up = 0;
		int down = 0;
		for (String gene : transitoryGenes.get(stage)) {
			double d = transitionDeltas.get(geneProbesetMap.get(gene).iterator().next()).get(stage);
			if (d > 0) { // down-regulating
				++down;
			} else {
				++up;
			}
		}
		return down;
	}
	
	public void getTransitoryGenes() {
		for (int i = 0; i < nTransitions; ++i) {
			if (nReplicas > 3) {
				transitoryGenesNCn429.put(i, new HashSet<String>());
			}
			else {
				transitoryGenesNC.put(i, new HashSet<String>());
				transitoryGenes429.put(i, new HashSet<String>());
			}
		}
		
		boolean once = false;
		ArrayList<Double> deltasNCn429 = null;
		ArrayList<Double> deltas429 = null;
		ArrayList<Double> deltasNC = null;
		for (String probeSetId: expData429.keySet()) {
			if (nReplicas > 3) {
				deltasNCn429 = transitionDeltasNCn429.get(probeSetId);
			}
			else {
				deltas429 = transitionDeltas429.get(probeSetId);
				deltasNC = transitionDeltasNC.get(probeSetId);
			}
			
			double wS[] = new double[]{0.5, 0.5, 1, 1.5, 2, 2.5, 3, 3.5};
			for (int stage = 0; stage < nStages - 1; ++stage) {
			
				double w = 3.5;
				double thresholdNC = transitionThreshold;
				double threshold429 = transitionThreshold;
				double thresholdNCn429 = transitionThreshold;
				
				if (nReplicas > 3) {
//					thresholdNCn429 = w * (stageNoiseStDNCn429[stage] + stageNoiseStDNCn429[stage + 1]);
					
					double sumV = stageNoiseStD429[stage] * stageNoiseStD429[stage] +
							      stageNoiseStD429[stage + 1] * stageNoiseStD429[stage + 1] + 
							      stageNoiseStDNC[stage] * stageNoiseStDNC[stage] + 
							      stageNoiseStDNC[stage + 1] * stageNoiseStDNC[stage + 1];
					if (stage == 0) {
						sumV = stageNoiseStD429[stage + 1] * stageNoiseStD429[stage + 1] 
								+ stageNoiseStDNC[stage + 1] * stageNoiseStDNC[stage + 1];
					}
					thresholdNCn429 = w * Math.sqrt(sumV);
					
					if (Math.abs(deltasNCn429.get(stage)) > thresholdNCn429) {
						transitoryGenesNCn429.get(stage).add(probesetGeneMap.get(probeSetId));
					}
				}
				else {
					threshold429 = w * (stageNoiseStD429[stage] + stageNoiseStD429[stage + 1]);
					if (Math.abs(deltas429.get(stage)) > threshold429) {
						transitoryGenes429.get(stage).add(probesetGeneMap.get(probeSetId));
					}
				
					thresholdNC = w * (stageNoiseStDNC[stage] + stageNoiseStDNC[stage + 1]);
					if (Math.abs(deltasNC.get(stage)) > thresholdNC) {
						transitoryGenesNC.get(stage).add(probesetGeneMap.get(probeSetId));
					}
				}
				
				if (once == true) {
					if (nReplicas > 3) {
						System.out.println(thresholdNCn429);
					}
					else {
						System.out.println(threshold429);
						System.out.println(thresholdNC);
					}
				}
			}
			once = false;
		}
		
		for (int i = 0; i < nTransitions; ++i) {
			if (nReplicas > 3) {
				double downNCn429 = getUpDownRegulation(i, transitoryGenesNCn429, transitionDeltasNCn429);
				System.out.print(
						(transitoryGenesNCn429.get(i).size() * 1.0 / geneProbesetMap.size()) 
						+ "\t" + transitoryGenesNCn429.get(i).size()); 
				compareEMTGenes(i, transitoryGenesNCn429, "Ncn429");
			}
			else {
				double down429 = getUpDownRegulation(i, transitoryGenes429, transitionDeltas429);
				double downNC = getUpDownRegulation(i, transitoryGenesNC, transitionDeltasNC);
				System.out.println(transitoryGenes429.get(i).size() + "\t" + (down429 / transitoryGenes429.get(i).size())
					+ "\t" + transitoryGenesNC.get(i).size() + "\t" + (downNC / transitoryGenesNC.get(i).size()));
			
//				if (transitoryGenes429.get(i).size() < 30) {
//					for (String s : transitoryGenes429.get(i)) {
//						System.out.println(s);
//					}
//					System.out.println("## ## ##");
//				}
//
//				if (transitoryGenesNC.get(i).size() < 65) {
//					int oe = 0;
//					for (String s : transitoryGenesNC.get(i)) {
//						if (oe % 2 == 0) {
//							System.out.print(s + " & ");
//						} else {
//							System.out.println(s + "\\\\");
//						}
//						++oe;
//					}
//					System.out.println("** ** **");
//				}
				
				compareEMTGenes(i, transitoryGenes429, "m429");
				compareEMTGenes(i, transitoryGenesNC, "NC");
			}
		}
		
		for (String s: emtGenes) {
//			System.out.print(s + "\t");
//			if (nReplicas > 3) {
//				int v = 0;
//				if (emtGeneTransitionNCn429.containsKey(s)) v = emtGeneTransitionNCn429.get(s);
//				System.out.println(v);
//			}
//			else {
//				int v = 0;
//				if (emtGeneTransition429.containsKey(s)) v = emtGeneTransition429.get(s);
//				System.out.print(v + "\t");
//				v = 0;
//				if (emtGeneTransitionNC.containsKey(s)) v = emtGeneTransitionNC.get(s);
//				System.out.println(v);
//			}
		}
	}
	
	private void updateEMTTransitionCount(String id, String gene) {
		HashMap<String, Integer> dummy = null;
		if (id.equals("NC")) dummy = emtGeneTransitionNC;
		if (id.equals("m429")) dummy = emtGeneTransition429;
		if (id.equals("Ncn429")) dummy = emtGeneTransitionNCn429;
		
		if (dummy.containsKey(gene)) {
			dummy.put(gene, 1 + dummy.get(gene));
		}
		else {
			dummy.put(gene, 1);
		}
	}
	
	private void compareEMTGenes(int transition, HashMap<Integer, HashSet<String>> transitoryGenes, String id) {
		double kount = 0;
		for (String gene: transitoryGenes.get(transition)) {
			if (emtGenes.contains(gene)) {
				++kount;
				updateEMTTransitionCount(id, gene);
				System.out.println(gene);
			}
		}
//		System.out.println("For " + id + " and transition " + transition + " contains emt " + (kount / emtGenes.size()));
		System.out.println(
				"\t" + (kount / emtGenes.size())
				+ "\t" + kount);
	}
	
	private void getGeneToProbesetHistogram() {
		TreeMap<Integer, Integer> geneToProbesetHist = new TreeMap();
		for (String g: geneProbesetMap.keySet()) {
			int sz = geneProbesetMap.get(g).size();
			if (geneToProbesetHist.containsKey(sz)) {
				geneToProbesetHist.put(sz, geneToProbesetHist.get(sz) + 1);
			}
			else {
				geneToProbesetHist.put(sz, 1);
			}
			
			if (geneProbesetMap.get(g).size() > 5000) { System.out.println(g); }
		}
		
		for (int i: geneToProbesetHist.keySet()) {
			System.out.println(i + "\t" + geneToProbesetHist.get(i));
		}
	}
	
	private void findGenes(String f) {
		for (String s: geneProbesetMap.keySet()) {
			if (s.contains(f)) {
				for (String r: geneProbesetMap.get(s)) {
					System.out.println(s + "\t" + r);
				}
			}
		}
	}
	
	public static void main(String[] args) throws Exception {
		CancerData cancerData = new CancerData();
		cancerData.loadCancerData();
		
//		cancerData.getGeneToProbesetHistogram();
//		cancerData.getCDFExpVal();
		
//		cancerData.normalizeData();
//		cancerData.getAmbientNoises();
		
		cancerData.nReplicas = 3;
		cancerData.getStageNoiseDistribution();

		cancerData.nReplicas = 9;
		cancerData.getSignificantStageValues(cancerData.ncN429);
//		cancerData.getTransitionDeltas();
//		cancerData.getCDFDeltas();
		
//		cancerData.checkProbesetConsistency();

//		cancerData.getTransitoryGenes();
		
//		cancerData.findGenes();
	}
}

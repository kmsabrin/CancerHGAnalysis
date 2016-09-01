import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;


public class CancerData {
	HashMap<String, ArrayList<Double>> geneExpressionValues429;
	HashMap<String, ArrayList<Double>> geneExpressionValuesNC;
	HashMap<String, HashSet<String>> geneProbesetMap;
	int[] hourStartIndices = {0, 3, 6, 9, 12, 15};
	int dataSize;
	int nReplicas = 3;
	int nTimePoints = 6;
	int nTransitions = 5;
	
	public CancerData() {
		geneExpressionValues429 = new HashMap();
		geneExpressionValuesNC = new HashMap();
		geneProbesetMap = new HashMap();
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
				insertMappedList(geneExpressionValues429, probeSetId, Double.parseDouble(tokens[i]));
				insertMappedList(geneExpressionValuesNC, probeSetId, Double.parseDouble(tokens[i]));
			}
			
			for (int i = 5; i <= 19; ++i) {
				insertMappedList(geneExpressionValues429, probeSetId, Double.parseDouble(tokens[i]));
			}
			
			for (int i = 20; i <= 34; ++i) {
				insertMappedList(geneExpressionValuesNC, probeSetId, Double.parseDouble(tokens[i]));
			}
			
			if (geneProbesetMap.containsKey(geneName)) {
				geneProbesetMap.get(geneName).add(probeSetId);
			}
			else {
				HashSet<String> hset = new HashSet();
				hset.add(probeSetId);
				geneProbesetMap.put(geneName, hset);
			}
			
			++dataSize;
		}
		
		scanner.close();
	}
	
	public void checkOffData(HashMap<String, ArrayList<Double>> expValues) {
		for (int hourIndex: hourStartIndices) {
			double leftGapValues[] = new double[dataSize];
			double rightGapValues[] = new double[dataSize];
			int idx1 = 0;
			for (String probeSetId: expValues.keySet()) {
				double replicatedExpValues[] = new double[nReplicas];
				int idx2 = 0;
				for(double d: expValues.get(probeSetId).subList(hourIndex, hourIndex + nReplicas)) {
					replicatedExpValues[idx2++] = d;
				}
				Arrays.sort(replicatedExpValues);
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
			double spearmanR429 = getSpearmanRhoTransitionDelta(sz, geneName, geneExpressionValues429);
			double spearmanRNC = getSpearmanRhoTransitionDelta(sz, geneName, geneExpressionValuesNC);
			System.out.println(spearmanR429 + "\t" + spearmanRNC);
		}
	}
	
	public void checkDataConsistency() {
		//checkOffData(geneExpressionValues429);
		//System.out.println("\n\n");
		//checkOffData(geneExpressionValuesNC);
		
		checkProbesetConsistency();
	}
	
	public static void main(String[] args) throws Exception {
		CancerData cancerData = new CancerData();
		cancerData.loadCancerData();
		cancerData.checkDataConsistency();
	}
}

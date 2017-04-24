import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.math3.stat.StatUtils;


public class TwoDataAnalysis {
	public static HashMap<String, Double> oldData = new HashMap();
	public static HashMap<String, Double> newData = new HashMap();
	public static HashMap<String, Double> diffData = new HashMap();
	
	public static void checkConsistency1() throws Exception {
		int oldIndex = 23; // nc=23, 429=8
		int newIndex = 4; // nc=4, 429=10
		Scanner scanner = new Scanner(new File("CancerData//5h_dataset.txt"));
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			String probeSetId = tokens[0];
			String gene = tokens[1];
			double avg = (Double.parseDouble(tokens[newIndex]) + Double.parseDouble(tokens[newIndex + 1])) / 2.0;
			newData.put(probeSetId, avg);
		}
		scanner.close();	
		
		scanner = new Scanner(new File("CancerData//entire_dataset.txt"));
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			String probeSetId = tokens[0];
			String gene = tokens[1];
			double avg = (Double.parseDouble(tokens[oldIndex]) 
					+ Double.parseDouble(tokens[oldIndex + 1]) + Double.parseDouble(tokens[oldIndex + 2])) / 3.0;
			oldData.put(probeSetId, avg);
		}
		scanner.close();
		
		int kount = 0;
		for (String s: newData.keySet()) {
			if (oldData.containsKey(s)) {
				++kount;
			}
		}
		
		double diffArray[] = new double[kount];
		kount = 0;
		for (String s: newData.keySet()) {
			if (oldData.containsKey(s)) {
//				double diff = Math.abs(oldData.get(s) - newData.get(s));
				double diff = (oldData.get(s) - newData.get(s));
				diffArray[kount++] = diff;
			}
		}
		
		System.out.println(StatUtils.mean(diffArray) + "\t" + Math.sqrt(StatUtils.variance(diffArray)));
	}
	
	public static void checkConsistency2() throws Exception {
		HashMap<String, Double> diffData = new HashMap();
		int ncIndex = 4; // nc=23, 429=8
		int index429 = 10; // nc=4, 429=10
//		Scanner scanner = new Scanner(new File("CancerData//entire_dataset.txt"));
		Scanner scanner = new Scanner(new File("CancerData//5h_dataset.txt"));
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			String probeSetId = tokens[0];
			String gene = tokens[1];
			double avgNC = (
					Double.parseDouble(tokens[ncIndex]) + 
					Double.parseDouble(tokens[ncIndex + 1]) 
					//+ Double.parseDouble(tokens[ncIndex + 2])
					) / 2.0;
			double avg429 = (
					Double.parseDouble(tokens[index429]) + 
					Double.parseDouble(tokens[index429 + 1]) 
					//+ Double.parseDouble(tokens[index429 + 2])
					) / 2.0;
			diffData.put(probeSetId, Math.abs(avgNC - avg429));
		}
		scanner.close();
		
		double diffArray[] = new double[diffData.size()];
		int kount = 0;
		for (String s: diffData.keySet()) {
			diffArray[kount++] = diffData.get(s);
		}
		
		System.out.println(StatUtils.mean(diffArray) + "\t" + Math.sqrt(StatUtils.variance(diffArray)));
	}
	
	public static void checkConsistency3() throws Exception {
		HashMap<String, Double[]> oldNC = new HashMap();
		HashMap<String, Double[]> old429 = new HashMap();		
		int indexNC = 23; // nc=23, 429=8
		int index429 = 8; // nc=4, 429=10
		Scanner scanner = new Scanner(new File("CancerData//entire_dataset.txt"));
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			String probeSetId = tokens[0];
			String gene = tokens[1];
			Double nc[] = new Double[3];
			nc[0] = Double.parseDouble(tokens[indexNC + 0]);
			nc[1] = Double.parseDouble(tokens[indexNC + 1]);
			nc[2] = Double.parseDouble(tokens[indexNC + 2]);
			oldNC.put(probeSetId, nc);
			Double m429[] = new Double[3];
			m429[0] = Double.parseDouble(tokens[index429 + 0]);
			m429[1] = Double.parseDouble(tokens[index429 + 1]);
			m429[2] = Double.parseDouble(tokens[index429 + 2]);
			old429.put(probeSetId, m429);
		}
		scanner.close();
		
		HashMap<String, Double[]> newNC = new HashMap();
		HashMap<String, Double[]> new429 = new HashMap();
		indexNC = 4; // nc=23, 429=8
		index429 = 10; // nc=4, 429=10
		scanner = new Scanner(new File("CancerData//5h_dataset.txt"));
		scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String tokens[] = line.split("\t");
			String probeSetId = tokens[0];
			String gene = tokens[1];
			Double nc[] = new Double[2];
			nc[0] = Double.parseDouble(tokens[indexNC + 0]);
			nc[1] = Double.parseDouble(tokens[indexNC + 1]);
			newNC.put(probeSetId, nc);
			Double m429[] = new Double[2];
			m429[0] = Double.parseDouble(tokens[index429 + 0]);
			m429[1] = Double.parseDouble(tokens[index429 + 1]);
			new429.put(probeSetId, m429);
		}
		scanner.close();
		
		double kount = 0;
		double valid = 0;
		ArrayList<Double> bigDiff = new ArrayList();
		for (String s: oldNC.keySet()) {
			Double range[] = getRange(oldNC.get(s), old429.get(s));
			double median = getMedian(newNC.get(s), new429.get(s));
			if (median >= range[0] && median <= range[1]) {
				++valid;
			}
			else {
				bigDiff.add(Math.min(Math.abs(median - range[0]), Math.abs(median - range[1])));
			}
			++kount;
		}
		
		System.out.println(valid / kount);
		
		double bigDiffArray[] = new double[bigDiff.size()];
		int idx = 0;
		for (double d: bigDiff) {
			bigDiffArray[idx++] = d;
		}
		System.out.println(StatUtils.mean(bigDiffArray) + "\t" + Math.sqrt(StatUtils.variance(bigDiffArray)));
	}
	
	public static double getMedian(Double dnc[], Double d429[]) {
		int idx = 2;
		double v[] = new double[idx * idx];
		int k = 0;
		for (int i = 0; i < idx; ++i) {
			for(int j = 0; j < idx; ++j) {
				double m = d429[i] - dnc[j];
				v[k++] = m;
			}
		}
		
		return StatUtils.percentile(v, 50);
	}
	
	public static Double[] getRange(Double dnc[], Double d429[]) {
		Double d[] = new Double[2];
		int idx = 3;
		double v[] = new double[idx * idx];
		int k = 0;
		for (int i = 0; i < idx; ++i) {
			for(int j = 0; j < idx; ++j) {
				double m = d429[i] - dnc[j];
				v[k++] = m;
			}
		}
		
		Arrays.sort(v);
		d[0] = v[0];
		d[1] = v[v.length - 1];
		return d;
	}
	
	public static void main(String[] args) throws Exception {
		checkConsistency3();
	}
}

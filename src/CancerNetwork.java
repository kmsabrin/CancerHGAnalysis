import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class CancerNetwork {
	HashSet<String> tfGenes;
	HashSet<String> nonTFGenes;
	HashMap<String, HashSet<String>> trrustNetwork;
	HashMap<String, HashSet<String>> nameAlias;
	HashMap<String, String> replaceAlias;
	HashSet<String> miR429Target;
	CancerData stageData;
	
	public CancerNetwork() {
		tfGenes = new HashSet();
		nonTFGenes = new HashSet();
		trrustNetwork = new HashMap();
		nameAlias = new HashMap();
		replaceAlias = new HashMap();
		miR429Target = new HashSet();
	}
	
	private void loadTRRUSTNetwork() throws Exception {
		Scanner scan = new Scanner(new File("CancerData/trrust_rawdata.txt"));
		
		while (scan.hasNext()) {
			String src = scan.next();
			String tgt = scan.next();
			String type = scan.next();
			String id = scan.next();
			
			if (trrustNetwork.containsKey(src)) {
				trrustNetwork.get(src).add(tgt);
			}
			else {
				HashSet<String> hset = new HashSet();
				hset.add(tgt);
				trrustNetwork.put(src, hset);
			}
			
			tfGenes.add(src);
			nonTFGenes.add(tgt);
		}
		
//		System.out.println(tfGenes.size());
//		System.out.println(nonTFGenes.size());
		scan.close();
	}
	
	private void loadmiR429Target() throws Exception {
		Scanner scan = new Scanner(new File("CancerData/mir429_target.txt"));
		while (scan.hasNext()) {
			String tgt = scan.next();
			miR429Target.add(tgt);
		}
		scan.close();
	}
	
	private void getNameAliasHelper(String s) {
		for (int i: stageData.significantStageGenes.keySet()) {
			for (String r: stageData.significantStageGenes.get(i)) {
				String ss = s.toLowerCase();
				String rr = r.toLowerCase();
				String rrTokens[] = rr.split("\\s+");
				for (String rrT: rrTokens) {
					if (rrT.equals(ss)) {
						if (rrTokens.length == 1) {
							break;
						}
						if (nameAlias.containsKey(s)) {
							nameAlias.get(s).add(r);
						}
						else {
							HashSet<String> hset = new HashSet();
							hset.add(r);
							nameAlias.put(s, hset);
						}
						break;
					}
				}
			}
		}
	}
	
	private void getNameAlias() {
		for (String s: tfGenes) {
			getNameAliasHelper(s);
		}
		
		for (String s: nonTFGenes) {
			getNameAliasHelper(s);
		}
		
		for (String s: miR429Target) {
			getNameAliasHelper(s);
		}
		
		for (String s: nameAlias.keySet()) {
			for (String r: nameAlias.get(s)) {
				replaceAlias.put(r, s);
			}
		}
		
		for (int i: stageData.significantStageGenes.keySet()) {
			HashSet<String> rmv = new HashSet();
			HashSet<String> add = new HashSet();
			for (String r: stageData.significantStageGenes.get(i)) {
				if (replaceAlias.containsKey(r)) {
					rmv.add(r);
					add.add(replaceAlias.get(r));
				}
			}
			stageData.significantStageGenes.get(i).removeAll(rmv);
			stageData.significantStageGenes.get(i).addAll(add);
		}
	}
	
	private void getNetwork() {
		HashMap<String, Integer> inKount = new HashMap();
		HashMap<String, Integer> outKount = new HashMap();
		
		for (int i = 1; i < stageData.nStages; ++i) {
			for (String s: stageData.significantStageGenes.get(i)) {
				String s_tmp = s + "_" + i; 
				
				if (miR429Target.contains(s)) {
					if (inKount.containsKey(s_tmp)) {
						inKount.put(s_tmp, inKount.get(s_tmp) + 1);
					}
					else {
						inKount.put(s_tmp,  1);
					}
				}
				
				if (!inKount.containsKey(s_tmp)) continue;
				
				for (int j = i + 1; j < stageData.nStages; ++j) {
					for (String r: stageData.significantStageGenes.get(j)) {
						String r_tmp = r + "_" + j;
						if (trrustNetwork.containsKey(s) && trrustNetwork.get(s).contains(r)) {							
							if (inKount.containsKey(r_tmp)) {
								inKount.put(r_tmp, inKount.get(r_tmp) + 1);
							}
							else {
								inKount.put(r_tmp,  1);
							}
							
							if (outKount.containsKey(s_tmp)) {
								outKount.put(s_tmp, outKount.get(s_tmp) + 1);
							}
							else {
								outKount.put(s_tmp, 1);
							}
							
//							if (r_tmp.equals("ATP1B1_4")) {
//								System.out.println(s_tmp + "\t" + r_tmp);
//							}
						}
					}
				}
				
			}
		}
		
		HashMap<Integer, HashSet<String>> visibleStageNodes = new HashMap();
		for (int i = 1; i <= 5; ++i) {
			HashSet<String> hset = new HashSet();
			visibleStageNodes.put(i, hset);
		}
	
//		System.out.println(inKount.get("ATP1B1_4"));
//		System.out.println(outKount.get("ATP1B1_4"));
//		System.out.println(outKount.containsKey("ATP1B1_4"));
		
		HashMap<String, String> nodeIdLabel = new HashMap();
		HashSet<String> printed = new HashSet();
		for (int i = 1; i < stageData.nStages; ++i) {
			int redBox = 0;
			int redBoxTF = 0;
			for (String s: stageData.significantStageGenes.get(i)) {
				String s_tmp = s + "_" + i;
				if (!inKount.containsKey(s_tmp)) continue;
				for (int j = i + 1; j < stageData.nStages; ++j) {
					int blueBox = 0;
					int blueBoxTF = 0;
					for (String r: stageData.significantStageGenes.get(j)) {
						String r_tmp = r + "_" + j;
						if (trrustNetwork.containsKey(s) && trrustNetwork.get(s).contains(r)) {
							if (!outKount.containsKey(r_tmp) && inKount.get(r_tmp) == 1) {
								blueBox++;
								if (tfGenes.contains(r)) {
									++blueBoxTF;
								}
							}
							else 
							{
								System.out.println(s_tmp + " -> " + r_tmp + ";");
								visibleStageNodes.get(i).add(s_tmp);
								visibleStageNodes.get(j).add(r_tmp);
								
								nodeIdLabel.put(s_tmp, s);
								nodeIdLabel.put(r_tmp, r);
							}
						}
					}
					if (blueBox > 0) {
						System.out.println(s_tmp + " -> " + s_tmp + "n" + blueBox + "_" + j + ";");
//						System.out.println(s_tmp + "n" + blueBox  + "_" + j + " [shape=box, style=filled, fillcolor=blue];");
						String blueBoxLabel = blueBox + "g (" + blueBoxTF + " TF)";
						System.out.println(s_tmp + "n" + blueBox  + "_" + j + " [shape=septagon, label=\"" + blueBoxLabel + "\"];");
						
						visibleStageNodes.get(j).add(s_tmp + "n" + blueBox + "_" + j);
					}
				}
				if (miR429Target.contains(s)) {
					if (!outKount.containsKey(s_tmp) && inKount.get(s_tmp) == 1) {
						redBox++;
						if (tfGenes.contains(s)) {
							redBoxTF++;
						}
					}
					else {
						String label = s;
						if (tfGenes.contains(s)) {
							System.out.println(s_tmp +  " [shape=ellipse, style=filled, fillcolor=orangered, label=\"" + label + "\"];");
						}
						else {
							System.out.println(s_tmp +  " [shape=box, style=filled, fillcolor=orangered, label=\"" + label + "\"];");
						}
						visibleStageNodes.get(i).add(s_tmp);
						printed.add(s_tmp);
					}
				}
			}
			String redBoxLabel = redBox + "g (" + redBoxTF + " TF)";
			System.out.println("miR429" +  "_" + i + "n" + redBox + 
					" [shape=septagon, style=filled, fillcolor=orangered, label=\"" + redBoxLabel + "\"];");
			visibleStageNodes.get(i).add("miR429" +  "_" + i + "n" + redBox);
		}
		
		for (String s_tmp: nodeIdLabel.keySet()) {
			if (printed.contains(s_tmp)) continue;
			String label = nodeIdLabel.get(s_tmp);
			if (tfGenes.contains(label)) {
				System.out.println(s_tmp +  " [shape=ellipse, style=filled, fillcolor=greenyellow, label=\"" + label + "\"];");
			}
			else {
				System.out.println(s_tmp +  " [shape=box, style=filled, fillcolor=greenyellow, label=\"" + label + "\"];");
			}

		}
		
		for (int i: visibleStageNodes.keySet()) {
			System.out.print("{");
			System.out.print("rank=same; ");
			for (String s: visibleStageNodes.get(i)) {
				System.out.print(s + "; ");
			}
			System.out.println("}");
		}
	}
	
	public static void main(String[] args) throws Exception {
		CancerNetwork cancerNet = new CancerNetwork();
		cancerNet.loadTRRUSTNetwork();
		cancerNet.loadmiR429Target();
		
		cancerNet.stageData = new CancerData();
		cancerNet.stageData.loadCancerData();
		cancerNet.stageData.nReplicas = 3;
		cancerNet.stageData.getStageNoiseDistribution();
		cancerNet.stageData.nReplicas = 9;
		cancerNet.stageData.getSignificantStageValues(cancerNet.stageData.ncN429);
		
//		for (int i: cancerNet.stageData.significantStageGenes.keySet()) {
//			System.out.println(i + "\t" + cancerNet.stageData.significantStageGenes.get(i).size());
//		}
		
		cancerNet.getNameAlias();
		
		cancerNet.getNetwork();
	}
}

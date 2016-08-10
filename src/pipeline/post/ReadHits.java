package pipeline.post;

import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

public class ReadHits {

	private ConcurrentHashMap<Integer, ConcurrentHashMap<Integer, Vector<Hit>>> hitMap = new ConcurrentHashMap<Integer, ConcurrentHashMap<Integer, Vector<Hit>>>();

	public void add(Hit h, int gi, int frame) {

		if (!hitMap.containsKey(gi))
			hitMap.put(gi, new ConcurrentHashMap<Integer, Vector<Hit>>());
		ConcurrentHashMap<Integer, Vector<Hit>> frameMap = hitMap.get(gi);

		if (!frameMap.containsKey(frame))
			frameMap.put(frame, new Vector<Hit>());
		Vector<Hit> hits = frameMap.get(frame);

		hits.add(h);

	}

	public ConcurrentHashMap<Integer, ConcurrentHashMap<Integer, Vector<Hit>>> getHitMap() {
		return hitMap;
	}

	public Vector<Hit> getAllHits() {
		Vector<Hit> allHits = new Vector<Hit>();
		for (Integer gi : hitMap.keySet()) {
			for (int frame : hitMap.get(gi).keySet()) {
				for (Hit h : hitMap.get(gi).get(frame)) {
					allHits.add(h);
				}
			}
		}
		return allHits;
	}

	public void print() {
		for (Integer gi : hitMap.keySet()) {
			System.out.println(">GI: " + gi);
			for (int frame : hitMap.get(gi).keySet()) {
				System.out.println("\t>>Frame: " + frame);
				for (Hit h : hitMap.get(gi).get(frame)) {
					h.print("\t");
				}
			}
		}
	}

	public void freeFrameHits(int gi, int frame) {
		for (Hit h : hitMap.get(gi).get(frame))
			h.freeMemory();
		hitMap.get(gi).remove(frame);
	}

	public void freeGiHits(int gi) {
		for (int frame : hitMap.get(gi).keySet())
			freeFrameHits(gi, frame);
		hitMap.remove(gi);
	}

}

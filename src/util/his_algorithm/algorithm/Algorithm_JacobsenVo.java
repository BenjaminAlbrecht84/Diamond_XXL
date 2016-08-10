package util.his_algorithm.algorithm;

import java.io.RandomAccessFile;
import java.util.Vector;

import pipeline.post.Hit;
import pipeline.post.HitRun_Rater;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;
import util.his_algorithm.splay_tree.SplayNode;
import util.his_algorithm.splay_tree.SplayTree;

public class Algorithm_JacobsenVo {

	public Vector<Hit> run(Hit[] sequence, HitRun_Rater scorer, Frame_Direction dir, RandomAccessFile rafSAM, RandomAccessFile rafDAA,
			String readID) {

		SplayTree L = new SplayTree();
		double max = 0;
		SplayNode lastElement = null;

		for (int i = 0; i < sequence.length; i++) {

			Hit e_i = sequence[i];
			int s_i = e_i.getId();

			SplayNode s = L.prevNode(s_i);
			double v = s != null ? s.getWeight() : 0;
			SplayNode t = s != null ? L.nextNode(s.getId()) : L.getMin();
			double w = t != null ? t.getWeight() : 0;

			double hicWeight = cmpWeight(s, e_i, sequence, scorer, dir, rafSAM, rafDAA, readID);
			while (t != null) {
				// if (v + e_i.getWeight() < w)
				// break;
				if (hicWeight < w)
					break;
				SplayNode tPrime = L.nextNode(t.getId());
				L.delete(t);
				t = tPrime;
				w = t != null ? t.getWeight() : 0;
			}

			if (t == null || s_i < t.getId()) {
				// SplayNode x = L.insert(e_i.getId(), v + e_i.getWeight());
				SplayNode x = L.insert(e_i.getId(), hicWeight);
				Object[] info = { i, s };
				x.setInfo(info);
				if (x.getWeight() > max) {
					max = x.getWeight();
					lastElement = x;
				}
			}

		}

		Vector<Hit> hic = new Vector<Hit>();
		if (lastElement != null) {
			hic.add(sequence[(int) lastElement.getInfo()[0]]);
			SplayNode v = (SplayNode) lastElement.getInfo()[1];
			while (v != null) {
				hic.add(0, sequence[(int) v.getInfo()[0]]);
				v = (SplayNode) v.getInfo()[1];
			}
		}

		return hic;

	}

	private double cmpWeight(SplayNode v, Hit e_i, Hit[] sequence, HitRun_Rater scorer, Frame_Direction dir, RandomAccessFile rafSAM,
			RandomAccessFile rafDAA, String readID) {

		Vector<Hit> hic = new Vector<Hit>();
		extractHitRun(v, hic, sequence);
		hic.add(e_i);

		// int[] res = scorer.run(hic, dir);
		// int score = res[0];
		// int length = res[1];
		//
		// double refLength = hic.get(0).getRef_length();
		// double runLength = length;
		// int coverage = (int) Math.floor((runLength / refLength) * 100.);
		//
		// return (int) Math.round(Math.log(Math.pow((double) score, (double)
		// coverage)));

		return (int) scorer.run(hic, dir, rafSAM, rafDAA, readID)[2];
	}

	private void extractHitRun(SplayNode v, Vector<Hit> run, Hit[] sequence) {
		if (v != null) {
			run.add(sequence[(int) v.getInfo()[0]]);
			SplayNode w = (SplayNode) v.getInfo()[1];
			while (w != null) {
				run.add(0, sequence[(int) w.getInfo()[0]]);
				w = (SplayNode) w.getInfo()[1];
			}
		}
	}

}

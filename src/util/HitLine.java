package util;

import java.util.Vector;

import pipeline.post.Hit;

public class HitLine {

	private final double m = 1. / 3.;

	private double b;
	private int step;

	public HitLine(int step, Vector<Hit> hits) {
		this.step = step;
		initLine(hits);
	}

	public void initLine(Vector<Hit> hits) {

		StringBuffer x = new StringBuffer("c(");
		StringBuffer y = new StringBuffer("c(");

		double sumQ = 0, sumR = 0;
		for (Hit h : hits) {

			if (!h.equals(hits.firstElement())) {
				x.append(",");
				y.append(",");
			}

			int q = h.getId() * step + h.getQuery_start();
			int r = h.getRef_start();
			x.append(q);
			y.append(r);

			sumQ += h.getId() * step + h.getQuery_start();
			sumR += h.getRef_start();

		}

		double meanQ = sumQ / (double) hits.size();
		double meanR = sumR / (double) hits.size();

		b = meanR - m * meanQ;

	}

	public double getDistance(Hit h) {

		double q = h.getId() * step + h.getQuery_start();
		double r = h.getRef_start();
		double d = Math.abs(m * q - r + b) / Math.sqrt(Math.pow(m, 2));

		return d;

	}

}

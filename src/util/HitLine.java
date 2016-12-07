package util;

import java.util.Vector;

import pipeline.post.Hit;
import pipeline.post.mode_one.Alignment_Generator_inParallel.Frame_Direction;

public class HitLine {

	private StringBuffer log;
	private double m;
	private int step;
	private double b;

	public HitLine(int step, Hit h1, Hit h2, Frame_Direction dir) {
		this.m = dir == Frame_Direction.Positiv ? 1. / 3. : -(1. / 3.);
		this.step = step;
		Vector<Hit> hits = new Vector<Hit>();
		hits.add(h1);
		hits.add(h2);
		initLine(hits);
	}

	public HitLine(int step, Vector<Hit> hits, Frame_Direction dir) {
		this.m = dir == Frame_Direction.Positiv ? 1. / 3. : -(1. / 3.);
		this.step = step;
		initLine(hits);
	}

	private void initLine(Vector<Hit> hits) {

		StringBuffer x = new StringBuffer("c(");
		StringBuffer y = new StringBuffer("c(");

		double sumQ = 0, sumR = 0;
		for (Hit h : hits) {

			if (!h.equals(hits.firstElement())) {
				x.append(",");
				y.append(",");
			}

			int q = Math.abs(h.getId()) * step + h.getQuery_start() - 1;
			int r = h.getRef_start();
			x.append(q);
			y.append(r);

			sumQ += q;
			sumR += r;

		}

		double meanQ = sumQ / (double) hits.size();
		double meanR = sumR / (double) hits.size();

		b = meanR - m * meanQ;

		// if (m < 0) {
		log = new StringBuffer();
		log = log.append("\n" + "x<-" + x + ")\n");
		log = log.append("y<-" + y + ")\n");
		log = log.append("plot(x,y)\n");
		log = log.append("abline(" + b + "," + m + ")\n");
		log = log.append("abline(lm(y ~ x))\n");
		log = log.append("lm(y ~ x)\n");
		log = log.append(b + " " + m + "\n");
		// }

	}

	public double getDistance(Hit h, Frame_Direction frameDirection) {

		int offset = Math.abs(h.getId()) * step;
		int qStart = h.getQuery_start() + offset - 1;

		double q = (double) qStart;
		double r = (double) h.getRef_start();
		double d = Math.abs(m * q - r + b) / Math.sqrt(Math.pow(m, 2));

		return d;

	}

	public StringBuffer getLog() {
		return log;
	}

}

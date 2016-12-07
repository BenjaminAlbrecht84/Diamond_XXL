package pipeline.init;

public class Read {

	private String id, seq, qual;

	public Read(String id, String seq, String qual) {
		this.id = id;
		this.seq =seq;
		this.qual = qual;
	}

	public String getId() {
		return id;
	}

	public String getSeq() {
		return seq;
	}

	public String getQual() {
		return qual;
	}

	public int getLength() {
		return seq.length();
	}

}

package gov.ornl.ccs;

public enum AdiosBufferAllocWhen {
    UNKNOWN (0),
    NOW (1), 
    LATER (2);

    private int code;

    AdiosBufferAllocWhen(int code) {
        this.code = code;
    }

    public int getCode() {
        return this.code;
    }
}

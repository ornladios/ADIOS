package gov.ornl.ccs;

public enum AdiosFlag {
    UNKNOWN (0),
    YES (1), 
    NO (2);

    private int code;

    AdiosFlag(int code) {
        this.code = code;
    }

    public int getCode() {
        return this.code;
    }
}


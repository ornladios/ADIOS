package gov.ornl.ccs;

public enum AdiosDatatype {
    UNKNOWN(-1)             /* (SIZE) */
    ,BYTE(0)                 /* (1) */
    ,SHORT(1)                /* (2) */
    ,INTEGER(2)              /* (4) */
    ,LONG(4)                 /* (8) */
    
    ,UNSIGNED_BYTE(50)       /* (1) */
    ,UNSIGNED_SHORT(51)      /* (2) */
    ,UNSIGNED_INTEGER(52)    /* (4) */
    ,UNSIGNED_LONG(54)       /* (8) */
    
    ,REAL(5)                 /* (4) */
    ,DOUBLE(6)               /* (8) */
    ,LONG_DOUBLE(7)          /* (16) */
    
    ,STRING(9)               /* (?) */
    ,COMPLEX(10)             /* (8) */
    ,DOUBLE_COMPLEX(11)      /* (16) */;
        
    private int code;

    AdiosDatatype(int code) {
        this.code = code;
    }

    public int getCode() {
        return this.code;
    }
}


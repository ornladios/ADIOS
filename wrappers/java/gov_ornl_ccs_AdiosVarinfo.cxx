#include "gov_ornl_ccs_AdiosVarinfo.h"

#include <adios_read.h>
#include <iostream>

/*
 * Class:     gov_ornl_ccs_AdiosVarinfo
 * Method:    adios_inq_var
 * Signature: (JJLjava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosVarinfo_adios_1inq_1var
(JNIEnv * env, jobject obj, jlong gp, jstring varname)
{
    const char *str;
    jboolean isCopy;
    jclass cls;
    jfieldID fid;

    ADIOS_VARINFO* vp;

    str = env->GetStringUTFChars(varname, &isCopy);
    if (str == NULL) return -1; /* OutOfMemoryError already thrown */

    vp = adios_inq_var((ADIOS_GROUP *)gp, str);
    env->ReleaseStringUTFChars(varname, str);
    
    cls = env->GetObjectClass(obj);
    if (cls == NULL) {
        return -1;
    }

    fid = env->GetFieldID(cls, "vp", "J");
    env->SetLongField(obj, fid, (jlong)vp);

    fid = env->GetFieldID(cls, "grpid", "I");
    env->SetIntField(obj, fid, vp->grpid);

    fid = env->GetFieldID(cls, "varid", "I");
    env->SetIntField(obj, fid, vp->varid);

    fid = env->GetFieldID(cls, "type", "I");
    env->SetIntField(obj, fid, vp->type);

    fid = env->GetFieldID(cls, "ndim", "I");
    env->SetIntField(obj, fid, vp->ndim);

    fid = env->GetFieldID(cls, "timedim", "I");
    env->SetIntField(obj, fid, vp->timedim);

    jlongArray dims = env->NewLongArray(vp->ndim);
    env->SetLongArrayRegion(dims, 0, vp->ndim, (jlong *) vp->dims);    

    fid = env->GetFieldID(cls, "dims", "[J");
    env->SetObjectField(obj, fid, dims);

    if (vp->ndim == 0) {
        int size = adios_type_size(vp->type, NULL);
        jbyteArray value = env->NewByteArray(size);
        env->SetByteArrayRegion(value, 0, size, (jbyte *) vp->value);

        fid = env->GetFieldID(cls, "value", "[B");
        env->SetObjectField(obj, fid, value);
    } 

    return 0;
}

/*
 * Class:     gov_ornl_ccs_AdiosVarinfo
 * Method:    adios_free_varinfo
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosVarinfo_adios_1free_1varinfo
(JNIEnv * env, jobject obj)
{
    jclass cls;
    jfieldID fid;
    jlong vp;

    cls = env->GetObjectClass(obj);
    fid = env->GetFieldID(cls, "vp", "J");
    vp = env->GetLongField(obj, fid);

    adios_free_varinfo((ADIOS_VARINFO*)vp);

    env->SetLongField(obj, fid, 0);

    return 0;
}

/*
 * Class:     gov_ornl_ccs_AdiosVarinfo
 * Method:    adios_read_var_byid
 * Signature: (JI[J[J)[D
 */
JNIEXPORT jdoubleArray JNICALL Java_gov_ornl_ccs_AdiosVarinfo_adios_1read_1var_1byid
(JNIEnv * env, jobject obj, jlong gp, jint varid, jlongArray start, jlongArray count)
{
    jint len = env->GetArrayLength(count);

    if (len != env->GetArrayLength(start))
        return NULL;

    jlong *startarr = env->GetLongArrayElements(start, NULL);
    jlong *countarr = env->GetLongArrayElements(count, NULL);

    jlong ncount = 1;
    for (jint i = 0; i < len; i++)
        ncount *= countarr[i];

    jdoubleArray result = env->NewDoubleArray(ncount);
    jdouble *data = env->GetDoubleArrayElements(result, NULL);

    int64_t nbytes = adios_read_var_byid((ADIOS_GROUP*)gp, (int) varid, (uint64_t *)startarr, (uint64_t *)countarr, (void *)data);

    env->ReleaseLongArrayElements(start, startarr, 0);
    env->ReleaseLongArrayElements(count, countarr, 0);
    env->ReleaseDoubleArrayElements(result, data, 0);

    return result;
}

#define READ_ARRAY(JTYPE, CTYPE)                                        \
    /* std::cout << __FUNCTION__ << "..." << std::endl; */              \
    int result;                                                         \
    jint len = env->GetArrayLength(count);                              \
    if (len != env->GetArrayLength(start))  return 0;                   \
    jlong *startarr = env->GetLongArrayElements(start, NULL);           \
    jlong *countarr = env->GetLongArrayElements(count, NULL);           \
    jlong ncount = 1;                                                   \
    for (jint i = 0; i < len; i++) ncount *= countarr[i];               \
    JTYPE *data = env->Get ##CTYPE ##ArrayElements(val, NULL);          \
    result = adios_read_var_byid((ADIOS_GROUP*)gp, (int) varid, (uint64_t *)startarr, (uint64_t *)countarr, (void *)data); \
    env->ReleaseLongArrayElements(start, startarr, 0);                  \
    env->ReleaseLongArrayElements(count, countarr, 0);                  \
    env->Release ##CTYPE ##ArrayElements(val, data, 0);                 \
    return result;

#define FUNC_READ_ARRAY(FNAME, JTYPE, CTYPE)                            \
    JNIEXPORT jint JNICALL FNAME                                        \
    (JNIEnv * env, jobject obj, jlong gp, jint varid, jlongArray start, jlongArray count, JTYPE ##Array val) \
    {                                                                   \
        READ_ARRAY(JTYPE, CTYPE);                                       \
    }

/*
 * Class:     gov_ornl_ccs_AdiosVarinfo
 * Method:    adios_read
 * Signature: (JI[J[J[B)I
 */

FUNC_READ_ARRAY (
    Java_gov_ornl_ccs_AdiosVarinfo_adios_1read__JI_3J_3J_3B,
    jbyte,
    Byte)

FUNC_READ_ARRAY (
    Java_gov_ornl_ccs_AdiosVarinfo_adios_1read__JI_3J_3J_3I,
    jint,
    Int)

FUNC_READ_ARRAY (
    Java_gov_ornl_ccs_AdiosVarinfo_adios_1read__JI_3J_3J_3J,
    jlong,
    Long)

FUNC_READ_ARRAY (
    Java_gov_ornl_ccs_AdiosVarinfo_adios_1read__JI_3J_3J_3F,
    jfloat,
    Float)

FUNC_READ_ARRAY (
    Java_gov_ornl_ccs_AdiosVarinfo_adios_1read__JI_3J_3J_3D,
    jdouble,
    Double)

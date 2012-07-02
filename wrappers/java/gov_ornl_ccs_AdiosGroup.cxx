#include "gov_ornl_ccs_AdiosGroup.h"

#include <adios_read.h>

/*
 * Class:     gov_ornl_ccs_AdiosGroup
 * Method:    adios_gopen
 * Signature: (JLjava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosGroup_adios_1gopen
(JNIEnv * env, jobject obj, jlong fp, jstring grpname)
{
    const char *str;
    jboolean isCopy;
    jclass cls;
    jfieldID fid;

    ADIOS_GROUP* gp;

    str = env->GetStringUTFChars(grpname, &isCopy);
    if (str == NULL) return -1; /* OutOfMemoryError already thrown */

    gp = adios_gopen((ADIOS_FILE *)fp, str);

    env->ReleaseStringUTFChars(grpname, str);
    
    cls = env->GetObjectClass(obj);
    if (cls == NULL) {
        return -1;
    }

    fid = env->GetFieldID(cls, "gp", "J");
    env->SetLongField(obj, fid, (jlong)gp);

    fid = env->GetFieldID(cls, "gh", "J");
    env->SetLongField(obj, fid, gp->gh);

    fid = env->GetFieldID(cls, "grpid", "I");
    env->SetIntField(obj, fid, gp->grpid);

    fid = env->GetFieldID(cls, "vars_count", "I");
    env->SetIntField(obj, fid, gp->vars_count);

    fid = env->GetFieldID(cls, "attrs_count", "I");
    env->SetIntField(obj, fid, gp->attrs_count);

    
    jclass stringClass = env->FindClass("java/lang/String");

    jobjectArray var_namelist = env->NewObjectArray(gp->vars_count, stringClass, NULL);

    for (int i = 0; i < gp->vars_count; i++)
    {
        jstring var_name = env->NewStringUTF(gp->var_namelist[i]);
        env->SetObjectArrayElement(var_namelist, i, var_name);
    }

    fid = env->GetFieldID(cls, "var_namelist", "[Ljava/lang/String;");    
    env->SetObjectField(obj, fid, var_namelist);

    jobjectArray attr_namelist = env->NewObjectArray(gp->attrs_count, stringClass, NULL);

    for (int i = 0; i < gp->attrs_count; i++)
    {
        jstring attr_name = env->NewStringUTF(gp->attr_namelist[i]);
        env->SetObjectArrayElement(attr_namelist, i, attr_name);
    }

    fid = env->GetFieldID(cls, "attr_namelist", "[Ljava/lang/String;");    
    env->SetObjectField(obj, fid, attr_namelist);

    return 0;
}

/*
 * Class:     gov_ornl_ccs_AdiosGroup
 * Method:    adios_gclose
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosGroup_adios_1gclose
(JNIEnv * env, jobject obj)
{
    jclass cls;
    jfieldID fid;
    jlong gp;
    jint result;

    cls = env->GetObjectClass(obj);
    fid = env->GetFieldID(cls, "gp", "J");
    gp = env->GetLongField(obj, fid);

    result = adios_gclose((ADIOS_GROUP*)gp);

    env->SetLongField(obj, fid, 0);

    return result;
}

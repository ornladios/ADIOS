#include "gov_ornl_ccs_AdiosFile.h"

#include <adios_read.h>

/*
 * Class:     gov_ornl_ccs_AdiosFile
 * Method:    adios_fopen
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosFile_adios_1fopen
(JNIEnv * env, jobject obj, jstring path, jlong comm)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;

    const char *str;
    int result;
    jboolean isCopy;
    jclass cls;
    jfieldID fid;

    ADIOS_FILE* fp;

    str = env->GetStringUTFChars(path, &isCopy);
    if (str == NULL) return -1; /* OutOfMemoryError already thrown */

    fp = adios_fopen(str, (MPI_Comm) comm);
    env->ReleaseStringUTFChars(path, str);

    cls = env->GetObjectClass(obj);
    if (cls == NULL) {
        return -1;
    }

    fid = env->GetFieldID(cls, "fp", "J");
    env->SetLongField(obj, fid, (jlong)fp);
    
    fid = env->GetFieldID(cls, "fh", "J");
    env->SetLongField(obj, fid, fp->fh);

    fid = env->GetFieldID(cls, "groups_count", "I");
    env->SetIntField(obj, fid, fp->groups_count);

    fid = env->GetFieldID(cls, "vars_count", "I");
    env->SetIntField(obj, fid, fp->vars_count);

    fid = env->GetFieldID(cls, "attrs_count", "I");
    env->SetIntField(obj, fid, fp->attrs_count);

    fid = env->GetFieldID(cls, "tidx_start", "I");
    env->SetIntField(obj, fid, fp->tidx_start);

    fid = env->GetFieldID(cls, "ntimesteps", "I");
    env->SetIntField(obj, fid, fp->ntimesteps);

    fid = env->GetFieldID(cls, "version", "I");
    env->SetIntField(obj, fid, fp->version);

    fid = env->GetFieldID(cls, "file_size", "J");
    env->SetIntField(obj, fid, fp->file_size);

    fid = env->GetFieldID(cls, "endianness", "I");
    env->SetIntField(obj, fid, fp->endianness);

    jclass stringClass = env->FindClass("java/lang/String");
    jobjectArray group_namelist = env->NewObjectArray(fp->groups_count, stringClass, NULL);

    for (int i = 0; i < fp->groups_count; i++)
    {
        jstring group_name = env->NewStringUTF(fp->group_namelist[i]);
        env->SetObjectArrayElement(group_namelist, i, group_name);
    }

    fid = env->GetFieldID(cls, "group_namelist", "[Ljava/lang/String;");    
    env->SetObjectField(obj, fid, group_namelist);

    return 0;
}

/*
 * Class:     gov_ornl_ccs_AdiosFile
 * Method:    adios_fclose
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_AdiosFile_adios_1fclose
(JNIEnv * env, jobject obj)
{
    jclass cls;
    jfieldID fid;
    jlong fp;
    jint result;

    cls = env->GetObjectClass(obj);
    fid = env->GetFieldID(cls, "fp", "J");
    fp = env->GetLongField(obj, fid);

    result = adios_fclose((ADIOS_FILE*)fp);

    fid = env->GetFieldID(cls, "fp", "J");
    env->SetLongField(obj, fid, 0);

    return result;
}

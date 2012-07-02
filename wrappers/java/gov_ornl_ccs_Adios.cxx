#include "gov_ornl_ccs_Adios.h"

#include <adios.h>
#include <iostream>
#include <string.h>

#define STR(A) #A
#define CONCAT2(A, B) A ## B
#define CONCAT3(A, B, C) CONCAT2(CONCAT2(A, B), C)
#define GET_ARRAY_ELEMENT(TYPE) CONCAT3(Get, TYPE, ArrayElements)
#define WRITE_ONE                                               \
    /* std::cout << __FUNCTION__ << "..." << std::endl;*/       \
    int result;                                                 \
    const char *str_var_name;                                   \
    str_var_name = env->GetStringUTFChars(var_name, NULL);      \
    if (str_var_name == NULL) return -1;                        \
    result = adios_write ((int64_t) fh, str_var_name, &val);    \
    env->ReleaseStringUTFChars(var_name, str_var_name);         \
    return result;

#define WRITE_ARRAY(JTYPE, CTYPE)                                       \
    /* std::cout << __FUNCTION__ << " ..." << std::endl;*/               \
    int result;                                                         \
    const char *str_var_name;                                           \
    JTYPE *valarr = env->Get ##CTYPE ##ArrayElements(val, NULL);        \
    str_var_name = env->GetStringUTFChars(var_name, NULL);              \
    if (str_var_name == NULL) return -1;                                \
    result = adios_write ((int64_t) fh, str_var_name, (void *) valarr); \
    env->ReleaseStringUTFChars(var_name, str_var_name);                 \
    env->Release ##CTYPE ##ArrayElements(val, valarr, 0);               \
    return result;

#define FUNC_WRITE_ONE(FNAME, JTYPE)                                    \
    JNIEXPORT jint JNICALL FNAME                                        \
    (JNIEnv * env, jclass cls, jlong fh, jstring var_name, JTYPE val)   \
    {                                                                   \
        WRITE_ONE;                                                      \
    }

#define FUNC_WRITE_ARRAY(FNAME, JTYPE, CTYPE)                           \
    JNIEXPORT jint JNICALL FNAME                                        \
    (JNIEnv * env, jclass cls, jlong fh, jstring var_name, JTYPE ##Array val) \
    {                                                                   \
        WRITE_ARRAY(JTYPE, CTYPE);                                      \
    }

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_init
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL 
Java_gov_ornl_ccs_Adios_adios_1init 
(JNIEnv *env, jclass cls, jstring xml_fname)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;
    const char *str;
    int result;
    jboolean isCopy;

    str = env->GetStringUTFChars(xml_fname, &isCopy);
    if (str == NULL) {
        return -1; /* OutOfMemoryError already thrown */
    }

    result = adios_init(str);
    //std::cout << "result = " << result << std::endl;
    
    env->ReleaseStringUTFChars(xml_fname, str);
    return result;
}

JNIEXPORT jint JNICALL 
Java_gov_ornl_ccs_Adios_adios_1init_1noxml
(JNIEnv *env, jclass cls)
{
    return adios_init_noxml();
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_open
 * Signature: (Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;J)J
 */
JNIEXPORT jlong JNICALL Java_gov_ornl_ccs_Adios_adios_1open
(JNIEnv * env, jclass cls, jstring group_name, jstring file_name, jstring mode, jlong comm)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;

    const char *str_group_name;
    const char *str_file_name;
    const char *str_mode;
    int result;
    static int64_t fd_p = 0;

    str_group_name = env->GetStringUTFChars(group_name, NULL);
    if (str_group_name == NULL) return -1;

    str_file_name = env->GetStringUTFChars(file_name, NULL);
    if (str_file_name == NULL) return -1;

    str_mode = env->GetStringUTFChars(mode, NULL);
    if (str_mode == NULL) return -1;

    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    //std::cout << "[IN] str_group_name = " << str_group_name << std::endl;
    //std::cout << "[IN] str_file_name = " << str_file_name << std::endl;
    //std::cout << "[IN] str_mode = " << str_mode << std::endl;
    //std::cout << "[IN] comm = " << (long) comm << std::endl;
    result = adios_open(&fd_p, str_group_name, str_file_name, str_mode, &comm);
    //std::cout << "[OUT] fd_p = " << fd_p << std::endl;
    //std::cout << "[OUT] result = " << result << std::endl;

    env->ReleaseStringUTFChars(file_name, str_file_name);
    env->ReleaseStringUTFChars(group_name, str_group_name);
    env->ReleaseStringUTFChars(mode, str_mode);

    /*
    uint64_t total_size;
    std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    adios_group_size (fd_p, 92, (uint64_t *) &total_size);
    std::cout << "[OUT] total size = " << total_size << std::endl;
    */

    return (jlong) fd_p;
}

JNIEXPORT jlong JNICALL Java_gov_ornl_ccs_Adios_adios_1group_1size
(JNIEnv * env, jclass cls, jlong fh, jlong group_size)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;

    uint64_t total_size = 0;
    int result;

    //std::cout << "[IN] fh = " << (int64_t) fh << std::endl;
    //std::cout << "[IN] group_size = " << (uint64_t) group_size << std::endl;
    result = adios_group_size ((int64_t) fh, (uint64_t) group_size, (uint64_t *) &total_size);
    //std::cout << "[OUT] fh = " << (int64_t) fh << std::endl;
    //std::cout << "[OUT] total size = " << total_size << std::endl;
    //std::cout << "[OUT] result = " << result << std::endl;

    return total_size;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_write
 * Signature: (JLjava/lang/String;B)I
 * Signature: (JLjava/lang/String;I)I
 * Signature: (JLjava/lang/String;J)I
 * Signature: (JLjava/lang/String;F)I
 * Signature: (JLjava/lang/String;D)I
 */
FUNC_WRITE_ONE (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2B, 
    jbyte)

FUNC_WRITE_ONE (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2I,
    jint)

FUNC_WRITE_ONE (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2J,
    jlong)

FUNC_WRITE_ONE (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2F,
    jfloat)

FUNC_WRITE_ONE (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2D,
    jdouble)

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_write
 * Signature: (JLjava/lang/String;[B)I
 * Signature: (JLjava/lang/String;[I)I
 * Signature: (JLjava/lang/String;[J)I
 * Signature: (JLjava/lang/String;[F)I
 * Signature: (JLjava/lang/String;[D)I
 */
FUNC_WRITE_ARRAY (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2_3B,
    jbyte,
    Byte)

FUNC_WRITE_ARRAY (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2_3I,
    jint,
    Int)

FUNC_WRITE_ARRAY (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2_3J,
    jlong,
    Long)

FUNC_WRITE_ARRAY (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2_3F,
    jfloat,
    Float)

FUNC_WRITE_ARRAY (
    Java_gov_ornl_ccs_Adios_adios_1write__JLjava_lang_String_2_3D,
    jdouble,
    Double)

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_close
 * Signature: (J)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1close
(JNIEnv * env, jclass cls, jlong fh)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;
    return adios_close ((int64_t) fh);
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_finalize
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1finalize
(JNIEnv * env, jclass cls, jint id)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;
    return adios_finalize (id);
}

#ifdef ADIOS_USE_MPI
/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_mpi_init
 * Signature: ([Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1mpi_1init
(JNIEnv * env, jclass cls, jobjectArray args)
{
    int argc = env->GetArrayLength(args);
    char **argv;
    int result; 
    
    for (int i = 0; i < argc; i++) 
    {
        jstring jstr = (jstring) env->GetObjectArrayElement(args, i);
        const char *chr = env->GetStringUTFChars(jstr, 0);

        //argv[i] = (char *) malloc(env->GetStringLength(jstr) + 1);
        argv[i] = new char[env->GetStringLength(jstr) + 1];
        strcpy(argv[i], chr);

        env->ReleaseStringUTFChars(jstr, chr);
    }

    result = MPI_Init(&argc, &argv);

    for (int i = 0; i < argc; i++) 
    {
        //free(argv[i]);
        delete[] argv[i]; 
    }

    return result;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_mpi_comm_rank
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1mpi_1comm_1rank
(JNIEnv * env, jclass cls, jlong comm)
{
    int rank;
    MPI_Comm_rank((MPI_Comm) comm, &rank);
    
    return rank;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_mpi_comm_size
 * Signature: (I)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1mpi_1comm_1size
(JNIEnv * env, jclass cls, jlong comm)
{
    int size;
    MPI_Comm_size((MPI_Comm) comm, &size);
    
    return size;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_mpi_finalize
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1mpi_1finalize
(JNIEnv * env, jclass cls)
{
    return MPI_Finalize();
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_mpi_comm_world
 * Signature: ()I
 */
JNIEXPORT jlong JNICALL Java_gov_ornl_ccs_Adios_adios_1mpi_1comm_1world
(JNIEnv * env, jclass cls)
{
    return (jlong) MPI_COMM_WORLD;
}
#endif

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_open_and_set_group_size
 * Signature: (Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;JJ)J
 */
JNIEXPORT jlong JNICALL Java_gov_ornl_ccs_Adios_adios_1open_1and_1set_1group_1size
(JNIEnv * env, jclass cls, jstring group_name, jstring file_name, jstring mode, jlong group_size, jlong comm)
{
    // Doesn't work, either
    //jlong fh = Java_gov_ornl_ccs_Adios_adios_1open(env, cls, group_name, file_name, mode, comm);
    //Java_gov_ornl_ccs_Adios_adios_1group_1size(env, cls, fh, group_size);

    //std::cout << __FUNCTION__ << "..." << std::endl;

    const char *str_group_name;
    const char *str_file_name;
    const char *str_mode;
    int result;
    static int64_t fd_p = 0;

    str_group_name = env->GetStringUTFChars(group_name, NULL);
    if (str_group_name == NULL) return -1;

    str_file_name = env->GetStringUTFChars(file_name, NULL);
    if (str_file_name == NULL) return -1;

    str_mode = env->GetStringUTFChars(mode, NULL);
    if (str_mode == NULL) return -1;

    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    //std::cout << "[IN] str_group_name = " << str_group_name << std::endl;
    //std::cout << "[IN] str_file_name = " << str_file_name << std::endl;
    //std::cout << "[IN] str_mode = " << str_mode << std::endl;
    //std::cout << "[IN] comm = " << (long) comm << std::endl;
    result = adios_open(&fd_p, str_group_name, str_file_name, str_mode, &comm);
    //std::cout << "[OUT] fd_p = " << fd_p << std::endl;
    //std::cout << "[OUT] result = " << result << std::endl;

    env->ReleaseStringUTFChars(file_name, str_file_name);
    env->ReleaseStringUTFChars(group_name, str_group_name);
    env->ReleaseStringUTFChars(mode, str_mode);

    uint64_t total_size;
    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    adios_group_size (fd_p, group_size, (uint64_t *) &total_size);
    //std::cout << "[OUT] total size = " << total_size << std::endl;

    return (jlong) fd_p;
}

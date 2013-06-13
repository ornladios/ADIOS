#include "gov_ornl_ccs_Adios.h"

#include <adios.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>

#define STR(A) #A
#define CONCAT2(A, B) A ## B
#define CONCAT3(A, B, C) CONCAT2(CONCAT2(A, B), C)
#define GET_ARRAY_ELEMENT(TYPE) CONCAT3(Get, TYPE, ArrayElements)
#define WRITE_ONE(JTYPE)                                        \
    /* std::cout << __FUNCTION__ << "..." << std::endl; */      \
    int result;                                                 \
    const char *str_var_name;                                   \
    str_var_name = env->GetStringUTFChars(var_name, NULL);      \
    if (str_var_name == NULL) return -1;                        \
    result = adios_write ((int64_t) fh, str_var_name, &val);    \
    env->ReleaseStringUTFChars(var_name, str_var_name);         \
    return result;

#define WRITE_ARRAY(JTYPE, CTYPE)                                       \
    /* std::cout << __FUNCTION__ << " ..." << std::endl;*/              \
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
        WRITE_ONE(JTYPE);                                               \
    }

#define FUNC_WRITE_ARRAY(FNAME, JTYPE, CTYPE)                           \
    JNIEXPORT jint JNICALL FNAME                                        \
    (JNIEnv * env, jclass cls, jlong fh, jstring var_name, JTYPE ##Array val) \
    {                                                                   \
        WRITE_ARRAY(JTYPE, CTYPE);                                      \
    }

#define STR_ALLOC(var) \
    const char *str_##var = env->GetStringUTFChars(var, NULL); \
    if (str_##var == NULL) return -1;

#define STR_ALLOC2(var) \
    const char *str_##var = env->GetStringUTFChars(var, NULL); \
    if (str_##var == NULL) goto end;

#define STR_FREE(var) \
    env->ReleaseStringUTFChars(var, str_##var);


/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_init
 * Signature: (Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL 
Java_gov_ornl_ccs_Adios_adios_1init 
(JNIEnv *env, jclass cls, jstring xml_fname, jlong comm)
{
    //std::cout << __FUNCTION__ << "..." << std::endl;
    int result;
    jboolean isCopy;

    STR_ALLOC(xml_fname);

    result = adios_init(str_xml_fname, (MPI_Comm) comm);
    //std::cout << "result = " << result << std::endl;
    
    STR_FREE(xml_fname);

    return result;
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

    int result;
    int64_t fd_p = 0;

    STR_ALLOC(group_name);
    STR_ALLOC(file_name);
    STR_ALLOC(mode);

    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    //std::cout << "[IN] str_group_name = " << str_group_name << std::endl;
    //std::cout << "[IN] str_file_name = " << str_file_name << std::endl;
    //std::cout << "[IN] str_mode = " << str_mode << std::endl;
    //std::cout << "[IN] comm = " << (long) comm << std::endl;
    result = adios_open(&fd_p, str_group_name, str_file_name, str_mode, (MPI_Comm) comm);
    //std::cout << "[OUT] fd_p = " << fd_p << std::endl;
    //std::cout << "[OUT] result = " << result << std::endl;

    STR_FREE(group_name);
    STR_FREE(file_name);
    STR_FREE(mode);

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

    //std::cout << "[IN] fh = " << fh << std::endl;
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

    int result;
    int64_t fd_p = 0;

    STR_ALLOC(group_name);
    STR_ALLOC(file_name);
    STR_ALLOC(mode);

    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    //std::cout << "[IN] str_group_name = " << str_group_name << std::endl;
    //std::cout << "[IN] str_file_name = " << str_file_name << std::endl;
    //std::cout << "[IN] str_mode = " << str_mode << std::endl;
    //std::cout << "[IN] comm = " << (long) comm << std::endl;
    result = adios_open(&fd_p, str_group_name, str_file_name, str_mode, (MPI_Comm) comm);
    //std::cout << "[OUT] fd_p = " << fd_p << std::endl;
    //std::cout << "[OUT] result = " << result << std::endl;

    STR_FREE(group_name);
    STR_FREE(file_name);
    STR_FREE(mode);

    uint64_t total_size;
    //std::cout << "[IN] fd_p = " << (int64_t) fd_p << std::endl;
    adios_group_size (fd_p, group_size, (uint64_t *) &total_size);
    //std::cout << "[OUT] total size = " << total_size << std::endl;

    return (jlong) fd_p;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_init_noxml
 * Signature: ()I
 */
JNIEXPORT jint JNICALL 
Java_gov_ornl_ccs_Adios_adios_1init_1noxml
(JNIEnv *env, jclass cls, jlong comm)
{
    return adios_init_noxml((MPI_Comm) comm);
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_allocate_buffer
 * Signature: (IJ)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1allocate_1buffer
(JNIEnv * env, jclass cls, jint when, jlong size)
{
    return adios_allocate_buffer((ADIOS_BUFFER_ALLOC_WHEN) when, size);
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_declare_group
 * Signature: (Ljava/lang/String;Ljava/lang/String;I)J
 */
JNIEXPORT jlong JNICALL Java_gov_ornl_ccs_Adios_adios_1declare_1group
(JNIEnv * env, jclass cls, jstring name, jstring time_index, jint stats)
{
    int64_t id_p;
    int result;

    STR_ALLOC(name);
    STR_ALLOC(time_index);

    result = adios_declare_group(&id_p, str_name, str_time_index, (ADIOS_FLAG) stats);

    STR_FREE(name);
    STR_FREE(time_index);

    return (jlong) id_p;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_define_var
 * Signature: (JLjava/lang/String;Ljava/lang/String;ILjava/lang/String;Ljava/lang/String;Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1define_1var
(JNIEnv * env, jclass cls, jlong group_id, jstring name, jstring path, jint type, jstring dimensions, jstring global_dimensions, jstring local_offsets)
{
    int result;

    STR_ALLOC(name);
    STR_ALLOC(path);
    STR_ALLOC(dimensions);
    STR_ALLOC(global_dimensions);
    STR_ALLOC(local_offsets);

    result = adios_define_var((int64_t) group_id, str_name, str_path, (ADIOS_DATATYPES) type, str_dimensions, str_global_dimensions, str_local_offsets);
    
    STR_FREE(name);
    STR_FREE(path);
    STR_FREE(dimensions);
    STR_FREE(global_dimensions);
    STR_FREE(local_offsets);

    return result;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_define_attribute
 * Signature: (JLjava/lang/String;Ljava/lang/String;ILjava/lang/String;Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1define_1attribute
(JNIEnv * env, jclass cls, jlong group_id, jstring name, jstring path, jint type, jstring value, jstring var)
{
    int result;
    
    STR_ALLOC(name);
    STR_ALLOC(path);
    STR_ALLOC(value);
    STR_ALLOC(var);
    
    result = adios_define_attribute((int64_t) group_id, str_name, str_path, (ADIOS_DATATYPES) type, str_value, str_var);

    STR_FREE(name);
    STR_FREE(path);
    STR_FREE(value);
    STR_FREE(var);

    return result;
}

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_select_method
 * Signature: (JLjava/lang/String;Ljava/lang/String;Ljava/lang/String;)I
 */
JNIEXPORT jint JNICALL Java_gov_ornl_ccs_Adios_adios_1select_1method
(JNIEnv * env, jclass cls, jlong group_id, jstring method, jstring parameters, jstring base_path)
{
    int result;
    
    STR_ALLOC(method);
    STR_ALLOC(parameters);
    STR_ALLOC(base_path);
    
    result = adios_select_method((int64_t) group_id, str_method, str_parameters, str_base_path);

    STR_FREE(method);
    STR_FREE(parameters);
    STR_FREE(base_path);
    
    return result;
}

#define READ_ONE(JTYPE)                                                 \
    /* std::cout << __FUNCTION__ << " ..." << std::endl; */             \
    int result;                                                         \
    JTYPE val;                                                          \
    const char *str_var_name;                                           \
    uint64_t read_size = sizeof(JTYPE);                                 \
    str_var_name = env->GetStringUTFChars(var_name, NULL);              \
    if (str_var_name == NULL) return -1;                                \
    result = adios_read ((int64_t) fh, str_var_name, &val, sizeof(JTYPE)); \
    env->ReleaseStringUTFChars(var_name, str_var_name);                 \
    return val;

#define FUNC_READ_ONE(FNAME, JTYPE)                                     \
    JNIEXPORT JTYPE JNICALL FNAME                                       \
    (JNIEnv * env, jclass cls, jlong fh, jstring var_name)              \
    {                                                                   \
        READ_ONE(JTYPE);                                                \
    }

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_read_byte_value
 * Signature: (JLjava/lang/String;)B
 * Signature: (JLjava/lang/String;)I
 * Signature: (JLjava/lang/String;)J
 * Signature: (JLjava/lang/String;)F
 * Signature: (JLjava/lang/String;)D
 */

FUNC_READ_ONE(
    Java_gov_ornl_ccs_Adios_adios_1read_1byte_1value,
    jbyte);

FUNC_READ_ONE(
    Java_gov_ornl_ccs_Adios_adios_1read_1int_1value,
    jint);

FUNC_READ_ONE(
    Java_gov_ornl_ccs_Adios_adios_1read_1long_1value,
    jlong);

FUNC_READ_ONE(
    Java_gov_ornl_ccs_Adios_adios_1read_1float_1value,
    jfloat);

FUNC_READ_ONE(
    Java_gov_ornl_ccs_Adios_adios_1read_1double_1value,
    jdouble);



#define READ_ARRAY(JTYPE, CTYPE)                                        \
    /* std::cout << __FUNCTION__ << " ..." << std::endl; */             \
    int result;                                                         \
    const char *str_var_name;                                           \
    JTYPE *valarr = env->Get ##CTYPE ##ArrayElements(val, NULL);        \
    uint64_t read_size = env->GetArrayLength(val) * sizeof(JTYPE);      \
    str_var_name = env->GetStringUTFChars(var_name, NULL);              \
    if (str_var_name == NULL) return -1;                                \
    result = adios_read ((int64_t) fh, str_var_name, (void *) valarr, read_size); \
    env->ReleaseStringUTFChars(var_name, str_var_name);                 \
    env->Release ##CTYPE ##ArrayElements(val, valarr, 0);               \
    return result;

#define FUNC_READ_ARRAY(FNAME, JTYPE, CTYPE)                            \
    JNIEXPORT jint JNICALL FNAME                                        \
    (JNIEnv * env, jclass cls, jlong fh, jstring var_name, JTYPE ##Array val) \
    {                                                                   \
        READ_ARRAY(JTYPE, CTYPE);                                      \
    }

/*
 * Class:     gov_ornl_ccs_Adios
 * Method:    adios_read
 * Signature: (JLjava/lang/String;[B)I
 * Signature: (JLjava/lang/String;[I)I
 * Signature: (JLjava/lang/String;[J)I
 * Signature: (JLjava/lang/String;[F)I
 * Signature: (JLjava/lang/String;[D)I
 */

FUNC_READ_ARRAY(
    Java_gov_ornl_ccs_Adios_adios_1read__JLjava_lang_String_2_3B,
    jbyte,
    Byte);

FUNC_READ_ARRAY(
    Java_gov_ornl_ccs_Adios_adios_1read__JLjava_lang_String_2_3I,
    jint,
    Int);

FUNC_READ_ARRAY(
    Java_gov_ornl_ccs_Adios_adios_1read__JLjava_lang_String_2_3J,
    jlong,
    Long);

FUNC_READ_ARRAY(
    Java_gov_ornl_ccs_Adios_adios_1read__JLjava_lang_String_2_3F,
    jfloat,
    Float);

FUNC_READ_ARRAY(
    Java_gov_ornl_ccs_Adios_adios_1read__JLjava_lang_String_2_3D,
    jdouble,
    Double);


